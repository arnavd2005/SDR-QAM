// qam_demo.cpp - command-line demo for qam_modem.h: text -> signal -> text.
//
// Build (macOS, Xcode command-line tools):
//   clang++ -std=c++17 -O2 qam_demo.cpp -o qam_demo
//
// Examples
//   ./qam_demo --M 16 --text "Hello QAM" --snr 25
//   ./qam_demo --M 64 --text "Hi" --clean --wav tx.wav     # save the signal
//   ./qam_demo --M 64 --rx tx.wav                          # decode a WAV
//   ./qam_demo --M 16 --text "Hello" --csv rx.csv          # symbols for plotting
//   ./qam_demo --sweep --text "The quick brown fox"        # BER table, all M
#include "qam_modem.h"

#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <random>

using namespace qam;

// ------------------------------------------------------------- channel
// Random delay (=> unknown timing and carrier phase), gain CH_GAIN, AWGN at Es/N0.
constexpr double CH_GAIN = 0.5;

static std::vector<double> channel(const std::vector<double>& s, double snr_db,
                                   size_t n_syms, std::mt19937& rng) {
    std::uniform_int_distribution<int> d(0, 3 * SPS - 1);
    const int delay = d(rng);
    std::vector<double> r(delay + s.size() + SPS, 0.0);
    double energy = 0;
    for (size_t i = 0; i < s.size(); ++i) {
        r[delay + i] = CH_GAIN * s[i];
        energy += r[delay + i] * r[delay + i];
    }
    const double es = energy / n_syms;                              // energy per symbol
    const double sigma = std::sqrt(es / (2 * std::pow(10.0, snr_db / 10)));
    std::normal_distribution<double> noise(0.0, sigma);
    for (double& v : r) v += noise(rng);
    return r;
}

// ------------------------------------------------------------- WAV I/O
static void put32(std::ofstream& f, uint32_t v) { f.write(reinterpret_cast<char*>(&v), 4); }
static void put16(std::ofstream& f, uint16_t v) { f.write(reinterpret_cast<char*>(&v), 2); }

static void write_wav(const std::string& path, const std::vector<double>& s) {
    std::ofstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot write " + path);
    const uint32_t data_bytes = uint32_t(s.size() * 2);
    f.write("RIFF", 4); put32(f, 36 + data_bytes); f.write("WAVE", 4);
    f.write("fmt ", 4); put32(f, 16); put16(f, 1); put16(f, 1);    // PCM, mono
    put32(f, FS); put32(f, FS * 2); put16(f, 2); put16(f, 16);
    f.write("data", 4); put32(f, data_bytes);
    for (double v : s) {
        const int16_t x = int16_t(std::lround(std::clamp(v, -1.0, 1.0) * 32767));
        f.write(reinterpret_cast<const char*>(&x), 2);
    }
}

static std::vector<double> read_wav(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot read " + path);
    std::vector<char> buf((std::istreambuf_iterator<char>(f)), {});
    auto u16 = [&](size_t i) { uint16_t v; std::memcpy(&v, &buf[i], 2); return v; };
    auto u32 = [&](size_t i) { uint32_t v; std::memcpy(&v, &buf[i], 4); return v; };
    if (buf.size() < 12 || std::memcmp(&buf[0], "RIFF", 4) || std::memcmp(&buf[8], "WAVE", 4))
        throw std::runtime_error("not a WAV file");
    int channels = 0, bits = 0;
    uint32_t rate = 0;
    for (size_t p = 12; p + 8 <= buf.size();) {           // walk chunks
        const uint32_t len = u32(p + 4);
        if (!std::memcmp(&buf[p], "fmt ", 4)) {
            channels = u16(p + 10); rate = u32(p + 12); bits = u16(p + 22);
        } else if (!std::memcmp(&buf[p], "data", 4)) {
            if (rate != FS || bits != 16 || channels < 1)
                throw std::runtime_error("need 16-bit PCM WAV at " + std::to_string(FS) + " Hz");
            const size_t frames = std::min<size_t>(len, buf.size() - p - 8) / (2 * channels);
            std::vector<double> x(frames);
            for (size_t i = 0; i < frames; ++i)
                x[i] = int16_t(u16(p + 8 + i * 2 * channels)) / 32768.0;   // first channel
            return x;
        }
        p += 8 + len + (len & 1);
    }
    throw std::runtime_error("no audio data in WAV");
}

// -------------------------------------------------------------- helpers
struct Run { RxResult rx; double ber; std::vector<double> tx; };

static Run run_once(int M, const std::string& text, bool clean, double snr, unsigned seed = 1) {
    Constellation c(M);
    std::mt19937 rng(seed);
    size_t nd = 0;
    Run out;
    out.tx = modulate(text, c, &nd);
    out.rx = demodulate(clean ? out.tx : channel(out.tx, snr, nd + PREAMBLE_LEN, rng), c);
    const auto ref = text_to_bits(text);
    if (out.rx.bits.size() >= ref.size()) {
        size_t errs = 0;
        for (size_t i = 0; i < ref.size(); ++i) errs += ref[i] != out.rx.bits[i];
        out.ber = double(errs) / ref.size();
    } else {
        out.ber = NAN;
    }
    return out;
}

static void usage() {
    std::printf(
        "usage: qam_demo [--M 4|8|16|32|64] [--text STR] [--snr dB | --clean]\n"
        "                [--wav out.wav] [--csv rx.csv] [--rx in.wav] [--sweep]\n"
        "  --snr    Es/N0 of the simulated channel (default 30)\n"
        "  --clean  no channel impairments\n"
        "  --wav    save the transmitted signal\n"
        "  --csv    save received symbols (i,q,ideal_i,ideal_q) for plotting\n"
        "  --rx     decode a 16-bit %g kHz WAV instead of simulating\n"
        "  --sweep  BER/CRC table for all M over a range of SNRs\n", FS / 1000.0);
}

// ----------------------------------------------------------------- main
int main(int argc, char** argv) {
    int M = 16;
    std::string text = "Hello from QAM!", wav, csv, rx;
    double snr = 30;
    bool clean = false, sweep = false;

    try {
        for (int i = 1; i < argc; ++i) {
            const std::string a = argv[i];
            auto next = [&]() -> std::string {
                if (i + 1 >= argc) { usage(); std::exit(1); }
                return argv[++i];
            };
            if (a == "--M") M = std::stoi(next());          // stoi/stod throw on bad numbers
            else if (a == "--text") text = next();
            else if (a == "--snr") snr = std::stod(next());
            else if (a == "--clean") clean = true;
            else if (a == "--wav") wav = next();
            else if (a == "--csv") csv = next();
            else if (a == "--rx") rx = next();
            else if (a == "--sweep") sweep = true;
            else { usage(); return a == "-h" || a == "--help" ? 0 : 1; }
        }

        Constellation c(M);                              // rejects an invalid M in every mode

        if (!rx.empty()) {
            const auto r = demodulate(read_wav(rx), c);
            std::printf("%d-QAM  CRC %s  MER %.1f dB  sync at sample %zu  phase corrected %+.1f deg\n"
                        "text: %s\n", M, r.frame.crc_ok ? "OK" : "FAIL", r.mer_db, r.start,
                        r.phase_deg, r.frame.text.c_str());
            return r.frame.crc_ok ? 0 : 2;
        }

        if (sweep) {
            const double snrs[] = {6, 9, 12, 15, 18, 21, 24, 27, 30};
            std::printf("text: '%s'\nBER (and CRC: + ok, x fail) vs Es/N0\n\n   M ", text.c_str());
            for (double s : snrs) std::printf("%9.0fdB", s);
            std::printf("\n");
            for (int m : {4, 8, 16, 32, 64}) {
                std::printf("  %2d ", m);
                for (double s : snrs) {
                    const auto r = run_once(m, text, false, s);
                    std::printf("  %8.4f%c", r.ber, r.rx.frame.crc_ok ? '+' : 'x');
                }
                std::printf("\n");
            }
            return 0;
        }

        const auto r = run_once(M, text, clean, snr);
        std::printf("%d-QAM | %d bits/symbol | %d baud = %d bit/s | signal %.0f ms\n", M, c.k,
                    BAUD, BAUD * c.k, 1000.0 * r.tx.size() / FS);
        if (clean) std::printf("channel: none\n");
        else std::printf("channel: Es/N0 %g dB, random delay+phase\n", snr);
        std::printf("sync at sample %zu, phase corrected %+.1f deg, MER %.1f dB, BER %.4f, CRC %s\n",
                    r.rx.start, r.rx.phase_deg, r.rx.mer_db, r.ber,
                    r.rx.frame.crc_ok ? "OK" : "FAIL");
        std::printf("sent:     %s\nreceived: %s\n", text.c_str(), r.rx.frame.text.c_str());

        if (!wav.empty()) { write_wav(wav, r.tx); std::printf("wrote %s\n", wav.c_str()); }
        if (!csv.empty()) {
            std::ofstream f(csv);
            f << "i,q,ideal_i,ideal_q\n";
            for (const auto& s : r.rx.symbols) {
                const auto p = c.points[c.slice(s)];
                f << s.real() << ',' << s.imag() << ',' << p.real() << ',' << p.imag() << '\n';
            }
            std::printf("wrote %s\n", csv.c_str());
        }
        return r.rx.frame.crc_ok ? 0 : 2;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
}
