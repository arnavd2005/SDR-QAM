#!/usr/bin/env python3
"""
qam_demo.py - generic M-QAM (4/8/16/32/64) modem demo: text -> signal -> text.

Signal chain
  TX: text -> frame [len(2B) | UTF-8 payload | CRC-16(2B)] -> bits
      -> Gray-mapped M-QAM symbols, prefixed by a known QPSK preamble
      -> upsample + root-raised-cosine (RRC) pulse shaping
      -> upconvert to a real audio-band carrier (default 6 kHz @ 48 kHz)
  CH: optional simulated channel: random delay (=> unknown timing AND carrier
      phase), gain, additive white Gaussian noise at a chosen Es/N0
  RX: downconvert -> RRC matched filter -> preamble correlation (timing)
      -> one-tap complex equalizer from preamble (gain + phase)
      -> nearest-point slicer -> bits -> parse frame -> check CRC -> text

Constellations (unit average energy, Gray-coded on each axis):
  4 = 2x2, 8 = 4x2, 16 = 4x4, 32 = 8x4, 64 = 8x8 (rectangular)

Only numpy (+ matplotlib for --plot). The real-valued audio-band signal and
table-based mapper were chosen so the same design ports to a microcontroller.

Examples
  python3 qam_demo.py --M 16 --text "Hello QAM" --snr 25 --plot qam16.png
  python3 qam_demo.py --M 64 --text "Hi" --wav tx.wav          # save signal
  python3 qam_demo.py --M 64 --rx tx.wav                       # decode a WAV
  python3 qam_demo.py --sweep --text "The quick brown fox"     # all M, many SNRs
"""
import argparse
import binascii
import sys
import wave

import numpy as np

# ---------------------------------------------------------------- parameters
FS = 48000            # sample rate (Hz)
BAUD = 2400           # symbol rate (symbols/s)
SPS = FS // BAUD      # samples per symbol (20)
FC = 6000             # carrier frequency (Hz)
BETA = 0.35           # RRC roll-off
SPAN = 8              # RRC length in symbols
PREAMBLE_LEN = 32     # known QPSK symbols for sync
GUARD = FS // 50      # 20 ms of silence before and after the signal
TX_PEAK = 0.9         # output peak, fraction of full scale
CH_GAIN = 0.5         # simulated channel gain


def rrc_taps(beta=BETA, sps=SPS, span=SPAN):
    """Root-raised-cosine FIR, normalized to unit energy."""
    t = np.arange(-span * sps // 2, span * sps // 2 + 1) / sps
    h = np.empty_like(t)
    for i, ti in enumerate(t):
        if abs(ti) < 1e-12:
            h[i] = 1 - beta + 4 * beta / np.pi
        elif abs(abs(ti) - 1 / (4 * beta)) < 1e-9:
            h[i] = beta / np.sqrt(2) * ((1 + 2 / np.pi) * np.sin(np.pi / (4 * beta))
                                       + (1 - 2 / np.pi) * np.cos(np.pi / (4 * beta)))
        else:
            h[i] = (np.sin(np.pi * ti * (1 - beta))
                    + 4 * beta * ti * np.cos(np.pi * ti * (1 + beta))) / \
                   (np.pi * ti * (1 - (4 * beta * ti) ** 2))
    return h / np.sqrt(np.sum(h ** 2))


RRC = rrc_taps()
PREAMBLE = np.exp(1j * (np.pi / 4 + np.pi / 2 *
                        np.random.default_rng(7).integers(0, 4, PREAMBLE_LEN)))


# ------------------------------------------------------------- constellation
def gray_axis(nbits):
    """Amplitude for each nbits value along one axis, Gray-coded: -(L-1)..(L-1)."""
    L = 1 << nbits
    amp = np.empty(L)
    for i in range(L):                    # i = position on the axis
        amp[i ^ (i >> 1)] = 2 * i - (L - 1)
    return amp


class QAM:
    def __init__(self, M):
        k = int(round(np.log2(M)))
        if 1 << k != M or not 2 <= k <= 6:
            raise ValueError("M must be 4, 8, 16, 32 or 64")
        self.M, self.k = M, k
        self.kI, self.kQ = (k + 1) // 2, k // 2      # bits on I / Q axes
        ai, aq = gray_axis(self.kI), gray_axis(self.kQ)
        mask = (1 << self.kQ) - 1
        pts = np.array([ai[s >> self.kQ] + 1j * aq[s & mask] for s in range(M)])
        self.points = pts / np.sqrt(np.mean(np.abs(pts) ** 2))   # Es = 1
        self.weights = 1 << np.arange(k - 1, -1, -1)

    def bits_to_symbols(self, bits):
        pad = (-len(bits)) % self.k
        bits = np.concatenate([bits, np.zeros(pad, np.uint8)])
        return self.points[bits.reshape(-1, self.k) @ self.weights]

    def slice(self, syms):
        """Nearest constellation index for each received symbol."""
        return np.argmin(np.abs(syms[:, None] - self.points[None, :]), axis=1)

    def indices_to_bits(self, idx):
        return ((idx[:, None] >> np.arange(self.k - 1, -1, -1)) & 1).astype(np.uint8).ravel()


# --------------------------------------------------------------------- frame
def text_to_bits(text):
    payload = text.encode("utf-8")
    if len(payload) > 65535:
        raise ValueError("text too long")
    crc = binascii.crc_hqx(payload, 0xFFFF)
    frame = len(payload).to_bytes(2, "big") + payload + crc.to_bytes(2, "big")
    return np.unpackbits(np.frombuffer(frame, np.uint8))


def bits_to_text(bits):
    """Returns (text, crc_ok, n_frame_bits) or (None, False, 0) if unparseable."""
    if len(bits) < 16:
        return None, False, 0
    n = int.from_bytes(np.packbits(bits[:16]).tobytes(), "big")
    need = 16 + 8 * n + 16
    if len(bits) < need:
        return None, False, 0
    raw = np.packbits(bits[:need]).tobytes()
    payload, crc = raw[2:2 + n], int.from_bytes(raw[2 + n:], "big")
    ok = binascii.crc_hqx(payload, 0xFFFF) == crc
    return payload.decode("utf-8", errors="replace"), ok, need


# --------------------------------------------------------------- modulator
def modulate(text, qam):
    data_syms = qam.bits_to_symbols(text_to_bits(text))
    syms = np.concatenate([PREAMBLE, data_syms])
    up = np.zeros(len(syms) * SPS, complex)
    up[::SPS] = syms
    bb = np.convolve(up, RRC)                                  # complex baseband
    n = np.arange(len(bb))
    s = np.sqrt(2) * np.real(bb * np.exp(2j * np.pi * FC * n / FS))  # real passband
    guard = np.zeros(GUARD)                                    # 20 ms silence
    s = np.concatenate([guard, s, guard])
    return TX_PEAK * s / np.max(np.abs(s)), data_syms


# ------------------------------------------------------------------ channel
def channel(s, snr_db, rng, n_syms):
    """Random delay (=> unknown timing and carrier phase), gain, AWGN at Es/N0."""
    delay = int(rng.integers(0, 3 * SPS))
    r = CH_GAIN * np.concatenate([np.zeros(delay), s, np.zeros(SPS)])
    es = np.sum(r ** 2) / n_syms                     # energy per symbol
    sigma = np.sqrt(es / (2 * 10 ** (snr_db / 10)))  # gives Es/N0 at the slicer
    return r + sigma * rng.standard_normal(len(r))


# -------------------------------------------------------------- demodulator
def demodulate(r, qam):
    n = np.arange(len(r))
    bb = np.sqrt(2) * r * np.exp(-2j * np.pi * FC * n / FS)    # downconvert
    y = np.convolve(bb, RRC)                                   # matched filter (also kills 2*FC)

    # timing: correlate against the preamble at symbol spacing
    tmpl = np.zeros((PREAMBLE_LEN - 1) * SPS + 1, complex)
    tmpl[::SPS] = PREAMBLE
    c = np.correlate(y, tmpl, "valid")                         # sum y[n+m]*conj(tmpl[m])
    start = int(np.argmax(np.abs(c)))
    h = c[start] / np.sum(np.abs(PREAMBLE) ** 2)               # complex gain+phase

    rx = y[start::SPS] / h                                     # equalized symbols
    data = rx[PREAMBLE_LEN:]
    idx = qam.slice(data)
    text, ok, nbits = bits_to_text(qam.indices_to_bits(idx))
    nsym = -(-nbits // qam.k) if nbits else len(data)
    data = data[:nsym]
    err = data - qam.points[idx[:nsym]]
    mer = 10 * np.log10(1 / max(np.mean(np.abs(err) ** 2), 1e-12))
    return dict(text=text, crc_ok=ok, symbols=data, all_symbols=rx[PREAMBLE_LEN:], mer_db=mer,
                start=start, phase_deg=np.degrees(np.angle(h)))


# ---------------------------------------------------------------- WAV I/O
def write_wav(path, s):
    with wave.open(path, "wb") as w:
        w.setnchannels(1); w.setsampwidth(2); w.setframerate(FS)
        w.writeframes((np.clip(s, -1, 1) * 32767).astype("<i2").tobytes())


def read_wav(path):
    with wave.open(path, "rb") as w:
        if w.getframerate() != FS or w.getsampwidth() != 2:
            raise ValueError(f"need 16-bit WAV at {FS} Hz")
        x = np.frombuffer(w.readframes(w.getnframes()), "<i2").astype(float) / 32768
        return x.reshape(-1, w.getnchannels())[:, 0]


# ------------------------------------------------------------------- plots
def plot(path, qam, tx, res, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(2, 2, figsize=(11, 8.5))
    fig.suptitle(title, fontsize=13)

    a = ax[0, 0]                                     # ideal constellation + bit labels
    a.scatter(qam.points.real, qam.points.imag, c="C0", s=40)
    for i, p in enumerate(qam.points):
        a.annotate(format(i, f"0{qam.k}b"), (p.real, p.imag), textcoords="offset points",
                   xytext=(0, 6), ha="center", fontsize=7 if qam.M >= 32 else 9)
    a.set_title(f"{qam.M}-QAM constellation (Gray)"); a.set_aspect("equal"); a.grid(alpha=.3)

    a = ax[0, 1]                                     # received constellation
    sy = res["symbols"]
    a.scatter(sy.real, sy.imag, s=8, alpha=.6, c="C1", label="received")
    a.scatter(qam.points.real, qam.points.imag, marker="x", c="k", s=40, label="ideal")
    a.set_title(f"Received symbols  (MER {res['mer_db']:.1f} dB)")
    a.set_aspect("equal"); a.grid(alpha=.3); a.legend(loc="upper right", fontsize=8)

    a = ax[1, 0]                                     # waveform
    t0 = GUARD
    seg = tx[t0:t0 + 12 * SPS]
    a.plot(np.arange(len(seg)) / FS * 1e3, seg, lw=1)
    a.set_title("Transmitted passband waveform (first 12 symbols)")
    a.set_xlabel("ms"); a.grid(alpha=.3)

    a = ax[1, 1]                                     # spectrum
    f = np.fft.rfftfreq(len(tx), 1 / FS)
    P = 20 * np.log10(np.abs(np.fft.rfft(tx * np.hanning(len(tx)))) + 1e-9)
    a.plot(f / 1e3, P - P.max(), lw=.8)
    a.set_xlim(0, FS / 2e3); a.set_ylim(-90, 5)
    a.set_title(f"Spectrum: carrier {FC/1e3:g} kHz, {BAUD} baud, roll-off {BETA}")
    a.set_xlabel("kHz"); a.set_ylabel("dB"); a.grid(alpha=.3)

    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


# -------------------------------------------------------------------- main
def run_once(M, text, snr, seed=1):
    qam = QAM(M)
    rng = np.random.default_rng(seed)
    tx, data_syms = modulate(text, qam)
    rx = tx if snr is None else channel(tx, snr, rng, len(data_syms) + PREAMBLE_LEN)
    res = demodulate(rx, qam)
    ref = text_to_bits(text)
    got = qam.indices_to_bits(qam.slice(res["all_symbols"]))[:len(ref)]
    res["ber"] = np.mean(got != ref[:len(got)]) if len(got) == len(ref) else float("nan")
    return qam, tx, res


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--M", type=int, default=16, choices=[4, 8, 16, 32, 64])
    ap.add_argument("--text", default="Hello from QAM!")
    ap.add_argument("--snr", type=float, default=30, help="Es/N0 in dB (use --clean for none)")
    ap.add_argument("--clean", action="store_true", help="no channel impairments")
    ap.add_argument("--wav", help="save transmitted signal to this WAV")
    ap.add_argument("--rx", help="decode this WAV instead of simulating")
    ap.add_argument("--plot", help="save a figure (PNG)")
    ap.add_argument("--sweep", action="store_true", help="test all M over a range of SNRs")
    a = ap.parse_args()

    if a.rx:
        qam = QAM(a.M)
        try:
            res = demodulate(read_wav(a.rx), qam)
        except (ValueError, OSError, wave.Error) as e:
            sys.exit(f"error: {e}")
        print(f"{a.M}-QAM  CRC {'OK' if res['crc_ok'] else 'FAIL'}  MER {res['mer_db']:.1f} dB  "
              f"sync at sample {res['start']}  phase corrected {res['phase_deg']:+.1f}°")
        print("text:", res["text"])
        return

    if a.sweep:
        snrs = [6, 9, 12, 15, 18, 21, 24, 27, 30]
        print(f"text: {a.text!r}\nBER (and CRC) vs Es/N0\n")
        print("   M  " + "".join(f"{s:>9}dB" for s in snrs))
        for M in (4, 8, 16, 32, 64):
            row = []
            for s in snrs:
                _, _, r = run_once(M, a.text, s)
                row.append(f"{r['ber']:9.4f}{'✓' if r['crc_ok'] else '✗'}")
            print(f"  {M:2d}  " + " ".join(row))
        return

    qam, tx, res = run_once(a.M, a.text, None if a.clean else a.snr)
    print(f"{a.M}-QAM | {qam.k} bits/symbol | {BAUD} baud = {BAUD*qam.k} bit/s | "
          f"signal {len(tx)/FS*1000:.0f} ms")
    print(f"channel: {'none' if a.clean else f'Es/N0 {a.snr:g} dB, random delay+phase'}")
    print(f"sync at sample {res['start']}, phase corrected {res['phase_deg']:+.1f}°, "
          f"MER {res['mer_db']:.1f} dB, BER {res['ber']:.4f}, CRC {'OK' if res['crc_ok'] else 'FAIL'}")
    print("sent:    ", a.text)
    print("received:", res["text"])
    if a.wav:
        write_wav(a.wav, tx); print("wrote", a.wav)
    if a.plot:
        plot(a.plot, qam, tx, res, f"{a.M}-QAM demo — Es/N0 "
             f"{'∞' if a.clean else f'{a.snr:g} dB'} — received: {res['text']!r}")
        print("wrote", a.plot)


if __name__ == "__main__":
    main()
