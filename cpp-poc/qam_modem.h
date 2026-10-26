// qam_modem.h - generic M-QAM (4/8/16/32/64) modem core, header-only C++17.
//
// Same design as qam_demo.py, and signal-compatible with it: a WAV made by one
// can be decoded by the other.
//
//   TX: text -> frame [len(2B) | UTF-8 payload | CRC-16(2B)] -> bits (MSB first)
//       -> Gray-mapped M-QAM symbols after a known 32-symbol QPSK preamble
//       -> root-raised-cosine pulse shaping -> real carrier (6 kHz @ 48 kHz)
//   RX: downconvert -> RRC matched filter -> preamble correlation (timing)
//       -> one-tap complex equalizer (gain + phase) -> nearest-point slicer
//       -> bits -> frame -> CRC check -> text
//
// No dependencies beyond the standard library. This is the desktop reference
// implementation; the Arduino port will reuse the same tables and logic in a
// streaming, fixed-point form.
#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace qam {

using cf = std::complex<double>;

// ------------------------------------------------------------- parameters
constexpr int    FS           = 48000;      // sample rate (Hz)
constexpr int    BAUD         = 2400;       // symbol rate (symbols/s)
constexpr int    SPS          = FS / BAUD;  // samples per symbol (20)
constexpr double FC           = 6000.0;     // carrier (Hz)
constexpr double BETA         = 0.35;       // RRC roll-off
constexpr int    SPAN         = 8;          // RRC length in symbols
constexpr int    PREAMBLE_LEN = 32;         // known QPSK symbols for sync
constexpr int    GUARD        = FS / 50;    // 20 ms of silence before and after the signal
constexpr double TX_PEAK      = 0.9;        // output peak, fraction of full scale
constexpr double PI           = 3.14159265358979323846;

// Preamble QPSK phases (index k -> angle pi/4 + k*pi/2); identical to the Python version.
static const uint8_t PREAMBLE_IDX[PREAMBLE_LEN] = {
    3, 2, 2, 3, 2, 3, 3, 0, 0, 1, 1, 3, 3, 0, 1, 3,
    0, 3, 0, 1, 3, 1, 1, 1, 2, 1, 3, 1, 1, 2, 2, 2};

// --------------------------------------------------------------- filters
inline std::vector<double> make_rrc(double beta = BETA, int sps = SPS, int span = SPAN) {
    const int half = span * sps / 2;
    std::vector<double> h(2 * half + 1);
    double energy = 0;
    for (int i = -half; i <= half; ++i) {
        const double t = double(i) / sps;
        double v;
        if (std::abs(t) < 1e-12) {
            v = 1 - beta + 4 * beta / PI;
        } else if (std::abs(std::abs(t) - 1 / (4 * beta)) < 1e-9) {
            v = beta / std::sqrt(2.0) * ((1 + 2 / PI) * std::sin(PI / (4 * beta)) +
                                         (1 - 2 / PI) * std::cos(PI / (4 * beta)));
        } else {
            v = (std::sin(PI * t * (1 - beta)) + 4 * beta * t * std::cos(PI * t * (1 + beta))) /
                (PI * t * (1 - (4 * beta * t) * (4 * beta * t)));
        }
        h[i + half] = v;
        energy += v * v;
    }
    for (double& v : h) v /= std::sqrt(energy);   // unit energy
    return h;
}

inline const std::vector<double>& rrc() {
    static const std::vector<double> h = make_rrc();
    return h;
}

inline const std::vector<cf>& preamble() {
    static const std::vector<cf> p = [] {
        std::vector<cf> v(PREAMBLE_LEN);
        for (int i = 0; i < PREAMBLE_LEN; ++i)
            v[i] = std::polar(1.0, PI / 4 + PI / 2 * PREAMBLE_IDX[i]);
        return v;
    }();
    return p;
}

// ------------------------------------------------------------------ CRC
// CRC-16/CCITT, poly 0x1021, init 0xFFFF (same as Python binascii.crc_hqx(data, 0xFFFF)).
inline uint16_t crc16(const uint8_t* d, size_t n, uint16_t crc = 0xFFFF) {
    for (size_t i = 0; i < n; ++i) {
        crc ^= uint16_t(d[i]) << 8;
        for (int b = 0; b < 8; ++b)
            crc = (crc & 0x8000) ? uint16_t((crc << 1) ^ 0x1021) : uint16_t(crc << 1);
    }
    return crc;
}

// -------------------------------------------------------- constellation
class Constellation {
public:
    int M, k, kI, kQ;              // order, bits/symbol, bits on I and Q axes
    std::vector<cf> points;        // index (bit pattern) -> unit-energy point

    explicit Constellation(int m) : M(m) {
        k = 0;
        while ((1 << k) < M) ++k;
        if ((1 << k) != M || k < 2 || k > 6)
            throw std::invalid_argument("M must be 4, 8, 16, 32 or 64");
        kI = (k + 1) / 2;
        kQ = k / 2;
        const auto ai = gray_axis(kI), aq = gray_axis(kQ);
        const int mask = (1 << kQ) - 1;
        points.resize(M);
        double es = 0;
        for (int s = 0; s < M; ++s) {
            points[s] = cf(ai[s >> kQ], aq[s & mask]);
            es += std::norm(points[s]);
        }
        es /= M;
        for (auto& p : points) p /= std::sqrt(es);
    }

    // bits (one per byte, values 0/1) -> symbols; zero-pads the last symbol
    std::vector<cf> map(const std::vector<uint8_t>& bits) const {
        std::vector<cf> out((bits.size() + k - 1) / k);
        for (size_t s = 0; s < out.size(); ++s) {
            int idx = 0;
            for (int b = 0; b < k; ++b) {
                const size_t j = s * k + b;
                idx = (idx << 1) | (j < bits.size() ? bits[j] : 0);
            }
            out[s] = points[idx];
        }
        return out;
    }

    // nearest constellation index
    int slice(cf y) const {
        int best = 0;
        double bd = std::norm(y - points[0]);
        for (int i = 1; i < M; ++i) {
            const double d = std::norm(y - points[i]);
            if (d < bd) { bd = d; best = i; }
        }
        return best;
    }

    void append_bits(int idx, std::vector<uint8_t>& bits) const {
        for (int b = k - 1; b >= 0; --b) bits.push_back(uint8_t((idx >> b) & 1));
    }

private:
    // amplitude for each nbits value along one axis, Gray-coded: -(L-1)..(L-1)
    static std::vector<double> gray_axis(int nbits) {
        const int L = 1 << nbits;
        std::vector<double> amp(L);
        for (int i = 0; i < L; ++i) amp[i ^ (i >> 1)] = 2 * i - (L - 1);
        return amp;
    }
};

// ---------------------------------------------------------------- frame
inline std::vector<uint8_t> text_to_bits(const std::string& text) {
    if (text.size() > 65535) throw std::invalid_argument("text too long");
    std::vector<uint8_t> frame;
    frame.push_back(uint8_t(text.size() >> 8));
    frame.push_back(uint8_t(text.size() & 0xFF));
    frame.insert(frame.end(), text.begin(), text.end());
    const uint16_t crc = crc16(reinterpret_cast<const uint8_t*>(text.data()), text.size());
    frame.push_back(uint8_t(crc >> 8));
    frame.push_back(uint8_t(crc & 0xFF));

    std::vector<uint8_t> bits;
    bits.reserve(frame.size() * 8);
    for (uint8_t byte : frame)
        for (int b = 7; b >= 0; --b) bits.push_back(uint8_t((byte >> b) & 1));
    return bits;
}

struct Frame {
    bool parsed = false;
    bool crc_ok = false;
    std::string text;
    size_t nbits = 0;              // bits the frame occupies
};

inline Frame bits_to_text(const std::vector<uint8_t>& bits) {
    Frame f;
    auto byte_at = [&](size_t i) {
        uint8_t v = 0;
        for (int b = 0; b < 8; ++b) v = uint8_t((v << 1) | bits[i * 8 + b]);
        return v;
    };
    if (bits.size() < 16) return f;
    const size_t n = (size_t(byte_at(0)) << 8) | byte_at(1);
    const size_t need = 16 + 8 * n + 16;
    if (bits.size() < need) return f;
    std::vector<uint8_t> payload(n);
    for (size_t i = 0; i < n; ++i) payload[i] = byte_at(2 + i);
    const uint16_t crc = uint16_t((byte_at(2 + n) << 8) | byte_at(3 + n));
    f.parsed = true;
    f.crc_ok = crc16(payload.data(), n) == crc;
    f.text.assign(payload.begin(), payload.end());
    f.nbits = need;
    return f;
}

// ------------------------------------------------------------ modulator
// Returns a real passband signal in [-TX_PEAK, TX_PEAK] with GUARD samples of silence on each side.
inline std::vector<double> modulate(const std::string& text, const Constellation& c,
                                    size_t* n_data_syms = nullptr) {
    const auto data = c.map(text_to_bits(text));
    if (n_data_syms) *n_data_syms = data.size();

    std::vector<cf> syms(preamble());
    syms.insert(syms.end(), data.begin(), data.end());

    const auto& h = rrc();
    std::vector<cf> bb(syms.size() * SPS + h.size() - 1);
    for (size_t m = 0; m < syms.size(); ++m)            // upsample + RRC
        for (size_t t = 0; t < h.size(); ++t) bb[m * SPS + t] += syms[m] * h[t];

    const size_t guard = GUARD;
    std::vector<double> s(bb.size() + 2 * guard, 0.0);
    double peak = 0;
    for (size_t n = 0; n < bb.size(); ++n) {             // upconvert
        const double ph = 2 * PI * FC * double(n) / FS;
        const double v = std::sqrt(2.0) * (bb[n].real() * std::cos(ph) - bb[n].imag() * std::sin(ph));
        s[guard + n] = v;
        peak = std::max(peak, std::abs(v));
    }
    for (double& v : s) v *= TX_PEAK / peak;
    return s;
}

// ---------------------------------------------------------- demodulator
struct RxResult {
    Frame frame;
    std::vector<cf> symbols;       // equalized data symbols in the frame
    std::vector<uint8_t> bits;     // all decided bits after the preamble
    double mer_db = 0;             // modulation error ratio
    size_t start = 0;              // sync position (samples)
    double phase_deg = 0;          // carrier phase that was corrected
};

inline RxResult demodulate(const std::vector<double>& r, const Constellation& c) {
    const auto& h = rrc();
    const auto& pre = preamble();

    // downconvert
    std::vector<cf> bb(r.size());
    for (size_t n = 0; n < r.size(); ++n)
        bb[n] = std::sqrt(2.0) * r[n] * std::polar(1.0, -2 * PI * FC * double(n) / FS);

    // matched filter (also removes the 2*FC mixing product)
    std::vector<cf> y(bb.size() + h.size() - 1);
    for (size_t n = 0; n < bb.size(); ++n)
        for (size_t t = 0; t < h.size(); ++t) y[n + t] += bb[n] * h[t];

    // timing: correlate with the preamble at symbol spacing
    const size_t span = size_t(PREAMBLE_LEN - 1) * SPS + 1;
    if (y.size() < span) throw std::runtime_error("signal too short");
    RxResult res;
    cf best_c = 0;
    double best = -1;
    for (size_t n = 0; n + span <= y.size(); ++n) {
        cf acc = 0;
        for (int m = 0; m < PREAMBLE_LEN; ++m) acc += y[n + size_t(m) * SPS] * std::conj(pre[m]);
        if (std::norm(acc) > best) { best = std::norm(acc); best_c = acc; res.start = n; }
    }
    const cf gain = best_c / double(PREAMBLE_LEN);       // |preamble|^2 = 1 per symbol
    res.phase_deg = std::arg(gain) * 180 / PI;

    // equalize + slice every symbol after the preamble
    std::vector<cf> eq;
    std::vector<int> idx;
    std::vector<uint8_t> bits;
    for (size_t n = res.start + size_t(PREAMBLE_LEN) * SPS; n < y.size(); n += SPS) {
        eq.push_back(y[n] / gain);
        idx.push_back(c.slice(eq.back()));
        c.append_bits(idx.back(), bits);
    }
    res.frame = bits_to_text(bits);

    const size_t nsym = res.frame.parsed ? (res.frame.nbits + c.k - 1) / c.k : eq.size();
    double err = 0;
    res.bits = bits;
    for (size_t i = 0; i < nsym; ++i) {
        res.symbols.push_back(eq[i]);
        err += std::norm(eq[i] - c.points[idx[i]]);
    }
    res.mer_db = 10 * std::log10(1.0 / std::max(err / std::max<size_t>(nsym, 1), 1e-12));
    return res;
}

}  // namespace qam
