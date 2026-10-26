# QAM Modem - Requirements and Design

```
Document   : README.md
Project    : Generic M-QAM text modem (Phase 1: desktop reference)
Version    : 1.0
Date       : 2026-10-25
Author     : Arnav Deshpande
Status     : Baseline for Phase 1; Phase 2 (Arduino) items marked [P2]
```

---

## 1. Purpose and scope

This document specifies a software modem that converts text into an
audio-band QAM signal and converts that signal back into text. It
supports 4, 8, 16, 32 and 64-QAM.

Phase 1 (this baseline) covers two desktop implementations that must
behave identically:

```
+----------------------+--------------------------------------------+
| Implementation       | Files                                      |
+----------------------+--------------------------------------------+
| Python reference     | py-poc/qam_demo.py                         |
|                      |   (numpy, optional matplotlib)             |
| C++17 implementation | cpp-poc/qam_modem.h                        |
|                      |   (modem core, header-only)                |
|                      | cpp-poc/qam_demo.cpp                       |
|                      |   (command-line harness)                   |
+----------------------+--------------------------------------------+
```

Phase 2 (out of scope here, listed in section 10.4) ports the modem core
to a microcontroller.

### 1.1 Requirement conventions

- Each requirement has a unique ID: `<AREA>-<NN>`.
- "shall" = mandatory; "should" = recommended.
- Verification methods: **T** = test, **A** = analysis, **I** = inspection.

```
Area codes
  SYS  system / parameters       FRM  framing
  MAP  constellation mapping     TX   transmitter
  CH   channel simulator         RX   receiver
  IO   files and CLI             PERF performance
  ARC  code architecture         P2   Phase 2 (future)
```

### 1.2 Running the Python demo

One-time setup, from the project root (matplotlib is only needed for
`--plot`):

```sh
python3 -m venv .venv
.venv/bin/pip install numpy matplotlib
```

Then run the demo with that Python (or `source .venv/bin/activate`
first and use plain `python`):

```sh
.venv/bin/python py-poc/qam_demo.py                       # 16-QAM, "Hello from QAM!", 30 dB
.venv/bin/python py-poc/qam_demo.py --M 64 --text "Hi" --snr 25
.venv/bin/python py-poc/qam_demo.py --M 16 --clean        # no noise or delay
.venv/bin/python py-poc/qam_demo.py --M 16 --snr 20 --plot qam16.png
.venv/bin/python py-poc/qam_demo.py --sweep --text "The quick brown fox jumps over the lazy dog"
```

Save a signal to a WAV and decode it again. The WAV files are
interchangeable with the C++ version (SYS-04). Decode with the same
`--M` that was used to send:

```sh
.venv/bin/python py-poc/qam_demo.py --M 64 --text "Hi" --wav tx.wav
.venv/bin/python py-poc/qam_demo.py --M 64 --rx tx.wav
```

All options are listed in 7.2 and by `--help`. To build the C++
version, see ARC-06.

---

## 2. System overview

The transmitter builds the signal in four stages. The channel distorts
it. The receiver undoes each stage in reverse order.

```
  TRANSMITTER
  +------------+    +------------+    +-------------+    +------------+
  | Frame text |--->| Map to QAM |--->| RRC shaping |--->| Upconvert  |
  | len,bytes, |    | Gray, k    |    | 20 samples/ |    | 6 kHz      |
  | CRC-16     |    | bits/symbol|    | symbol      |    | carrier    |
  +------------+    +------------+    +-------------+    +-----+------+
                                                               |
                                            real audio signal  |
                                            (WAV or simulated) v
                                                         +------------+
                                                         |  Channel   |
                                                         | delay,gain,|
                                                         | noise      |
                                                         +-----+------+
                                                               |
  RECEIVER                                                     v
  +------------+    +------------+    +-------------+    +------------+
  | Slice,     |<---|   Sync     |<---| Matched     |<---| Downconvert|
  | decode,    |    | timing,    |    | filter      |    | back to    |
  | check CRC  |    | gain,phase |    | (same RRC)  |    | I and Q    |
  +-----+------+    +------------+    +-------------+    +------------+
        |
        v
      text + diagnostics (MER, BER, sync position, phase)
```

### 2.1 Fixed system parameters

```
+--------------+---------+-------------------------------------------+
| Name         | Value   | Meaning                                   |
+--------------+---------+-------------------------------------------+
| FS           | 48000   | sample rate, Hz                           |
| BAUD         | 2400    | symbol rate, symbols/s                    |
| SPS          | 20      | samples per symbol (FS / BAUD)            |
| FC           | 6000    | carrier frequency, Hz                     |
| BETA         | 0.35    | root-raised-cosine roll-off               |
| SPAN         | 8       | RRC length in symbols (161 taps)          |
| PREAMBLE_LEN | 32      | known QPSK symbols used for sync          |
| M            | 4..64   | constellation order: 4, 8, 16, 32 or 64   |
+--------------+---------+-------------------------------------------+
```

| ID     | Requirement | Verif. |
|--------|-------------|--------|
| SYS-01 | The modem shall use the parameter values in 2.1. | I |
| SYS-02 | All parameters shall be defined in exactly one place per implementation. | I |
| SYS-03 | The modem shall support M = 4, 8, 16, 32 and 64 and reject any other value with an error. | T |
| SYS-04 | The Python and C++ implementations shall produce signals that each can decode from the other (WAV interoperability). | T |

---

## 3. Background: DSP concepts used

This section explains the design; it contains no requirements.

**Sampling.** The signal is stored as 48,000 values per second. Nyquist
limits representable frequencies to 24 kHz; the signal sits near 6 kHz.

**Symbols, not bits.** Each symbol carries k = log2(M) bits. At 2400
symbols/s:

```
  M      k (bits/symbol)   bit rate
  -----  ---------------   ----------
   4          2             4800 bit/s
   8          3             7200 bit/s
  16          4             9600 bit/s
  32          5            12000 bit/s
  64          6            14400 bit/s
```

**Complex numbers and I/Q.** Each symbol is a point I + jQ. A cosine and
a sine at the same frequency are orthogonal, so I rides on a cosine, Q
on a sine, and one real signal carries both.

**Filtering is convolution.** A filter is a list of taps; each output
sample is a weighted sum of nearby input samples.

**Pulse shaping and intersymbol interference (ISI).** A raised-cosine
pulse overlaps its neighbours but is zero at every other symbol's centre,
so sampling at the right instants sees no ISI:

```
  amplitude
     ^          symbol n
     |            _
     |          /   \
     |         /     \
   0 +--.--.--/---+---\--.--.--.---> time
     |   \_/ |    |    | \_/
     |       n-1  n   n+1
     |   pulse is 0 at n-1, n+1, n+2 ... -> no ISI at sample instants
```

**Matched filtering.** The raised cosine is split into two root-raised-
cosine (RRC) halves, one at each end. Together they form a raised
cosine, and the receiver half maximises SNR at the sampling instant.

**Noise and SNR.** Es/N0 is symbol energy over noise density, in dB.
Every 10 dB is a factor of 10.

**Correlation.** Sliding a known pattern along a signal produces a peak
where the pattern occurs. The receiver uses this to find the preamble.

---

## 4. Framing requirements

### 4.1 Frame format

```
  byte:  0       1       2 ............. 2+N-1   2+N     2+N+1
        +-------+-------+-------------------------+-------+-------+
        | LEN_H | LEN_L |   payload (UTF-8, N B)  | CRC_H | CRC_L |
        +-------+-------+-------------------------+-------+-------+
        |<- length N -->|                         |<- CRC-16 ---->|
            big-endian                               of payload only

  Total frame = N + 4 bytes = 8 * (N + 4) bits, transmitted MSB first.
  Example: "Hi" -> N = 2 -> 6 bytes -> 48 bits.
```

| ID     | Requirement | Verif. |
|--------|-------------|--------|
| FRM-01 | The payload shall be the text encoded as UTF-8. | T |
| FRM-02 | The frame shall be [2-byte big-endian length][payload][2-byte CRC], as in 4.1. | I, T |
| FRM-03 | The CRC shall be CRC-16/CCITT: poly 0x1021, init 0xFFFF, no reflection, no final XOR (identical to Python `binascii.crc_hqx(data, 0xFFFF)`). | T |
| FRM-04 | Payloads longer than 65535 bytes shall be rejected with an error. | T |
| FRM-05 | Bytes shall be serialised to bits most significant bit first. | T |
| FRM-06 | The receiver shall report whether the CRC matched, and shall still return the decoded text when it does not. | T |

---

## 5. Constellation mapping requirements

### 5.1 Grid shape

Bits per symbol are split between the axes: kI = ceil(k/2) on I and
kQ = floor(k/2) on Q.

```
  M     k   kI  kQ   grid (I x Q)   points per axis
  ---   -   --  --   ------------   ---------------
   4    2   1   1        2 x 2       +-1
   8    3   2   1        4 x 2       I: +-1,+-3   Q: +-1
  16    4   2   2        4 x 4       +-1, +-3
  32    5   3   2        8 x 4       I: +-1..+-7  Q: +-1,+-3
  64    6   3   3        8 x 8       +-1, +-3, +-5, +-7
                         (raw levels, before energy normalisation)
```

### 5.2 Example: 16-QAM with Gray labels

The symbol index is [I bits][Q bits]. Neighbouring points differ in
exactly one bit.

```
                        Q
                        ^
                        |
     0010      0110     |     1110      1010      Q = +0.949  (Q bits 10)
                        |
     0011      0111     |     1111      1011      Q = +0.316  (Q bits 11)
  ----------------------+-------------------------> I
     0001      0101     |     1101      1001      Q = -0.316  (Q bits 01)
                        |
     0000      0100     |     1100      1000      Q = -0.949  (Q bits 00)
                        |
   -0.949    -0.316        +0.316    +0.949       I level
   (I 00)    (I 01)        (I 11)    (I 10)       I bits
```

| ID     | Requirement | Verif. |
|--------|-------------|--------|
| MAP-01 | Constellations shall be rectangular grids per 5.1. | I, T |
| MAP-02 | Each axis shall be Gray-coded: axis position i maps to bit value i XOR (i >> 1). | T |
| MAP-03 | The symbol index shall be formed as (I bits << kQ) OR (Q bits). | I |
| MAP-04 | Points shall be scaled to unit average energy (mean of abs(point)^2 = 1) for every M. | T |
| MAP-05 | If the bit count is not a multiple of k, the last symbol shall be zero-padded. | T |
| MAP-06 | The receiver shall decide each symbol by minimum Euclidean distance (nearest point). | I |

### 5.3 Why larger M needs more SNR

```
  M     minimum distance between points (unit energy)
  ---   --------------------------------------------
   4    ~1.41   ##############
  16    ~0.63   ######
  64    ~0.31   ###
        Each extra pair of bits costs about 6 dB of SNR.
```

---

## 6. Transmitter, channel and receiver requirements

### 6.1 Transmitted signal layout

```
  |<- 20 ms ->|<-- 32 symbols -->|<----- N data symbols ----->|<-tail->|<- 20 ms ->|
  +-----------+------------------+----------------------------+--------+-----------+
  |  silence  | QPSK preamble    | QAM data (frame from 4.1)  | RRC    |  silence  |
  |  (zeros)  | (fixed pattern)  |                            | tail   |  (zeros)  |
  +-----------+------------------+----------------------------+--------+-----------+
               one symbol = 20 samples = 0.4167 ms
```

### 6.2 Occupied spectrum

```
  level
    ^
    |              +-----------+
    |             /|           |\        bandwidth = BAUD * (1 + BETA)
    |            / |           | \                 = 2400 * 1.35
    |           /  |           |  \                = 3240 Hz
  --+----------+---+-----+-----+---+-----------------------------> Hz
    0        4380        6000        7620                   24000
                        carrier
```

### 6.3 Preamble

The preamble is 32 QPSK symbols at angle pi/4 + k * pi/2, where k is
taken in order from this fixed table (shared by both implementations):

```
  3 2 2 3 2 3 3 0   0 1 1 3 3 0 1 3   0 3 0 1 3 1 1 1   2 1 3 1 1 2 2 2
```

### 6.4 Transmitter

| ID    | Requirement | Verif. |
|-------|-------------|--------|
| TX-01 | The transmitter shall prepend the preamble in 6.3 to the data symbols. | I, T |
| TX-02 | Symbols shall be upsampled by SPS (one symbol, then SPS-1 zeros) and convolved with the RRC filter. | I |
| TX-03 | The RRC filter shall have 161 taps (SPAN * SPS + 1), roll-off BETA, and unit energy. | T |
| TX-04 | The passband signal shall be s(n) = sqrt(2) * [ I(n) cos(2 pi FC n / FS) - Q(n) sin(2 pi FC n / FS) ]. | I |
| TX-05 | The output shall be scaled to a peak of 0.9 full scale. | T |
| TX-06 | The output shall have 20 ms of silence before and after the signal (6.1). | T |

### 6.5 Channel simulator (test only)

| ID    | Requirement | Verif. |
|-------|-------------|--------|
| CH-01 | The simulator shall add a random delay of 0 to 3*SPS - 1 samples (unknown timing and carrier phase). | I |
| CH-02 | The simulator shall scale the signal by a gain of 0.5. | I |
| CH-03 | The simulator shall add white Gaussian noise sized to the requested Es/N0 at the slicer. | A |
| CH-04 | The simulator shall be seedable so that runs are repeatable. | T |
| CH-05 | A "clean" mode shall pass the signal through with no impairments. | T |

Note: the carrier advances 45 degrees per sample (6000/48000 cycle), so
any delay also rotates the constellation by an unknown angle. The
receiver must correct this (RX-04).

### 6.6 Receiver

```
  Sync by correlation (RX-03)

  correlation
  magnitude
    ^                      peak = start of preamble
    |                         |
    |                         #
    |                         #
    |    .   .  .    .   .    #   .  .    .   .
    +--------------------------+------------------> sample offset
                            start
    complex value at the peak = channel gain + phase (used by RX-04)
```

| ID    | Requirement | Verif. |
|-------|-------------|--------|
| RX-01 | The receiver shall downconvert by multiplying by sqrt(2) * exp(-j 2 pi FC n / FS). | I |
| RX-02 | The receiver shall apply the same RRC filter as a matched filter (also removing the 2*FC product). | I |
| RX-03 | The receiver shall locate the frame by correlating against the preamble at symbol spacing and taking the peak magnitude. | T |
| RX-04 | The receiver shall correct gain and phase by dividing all samples by (peak value / PREAMBLE_LEN) (one-tap equalizer). | T |
| RX-05 | After the preamble, the receiver shall take one sample per symbol and slice it (MAP-06). | I |
| RX-06 | The receiver shall read the length field, then exactly the bits the frame needs, then check the CRC. | T |
| RX-07 | The receiver shall report: text, CRC status, equalized symbols, MER (dB), sync position (samples) and corrected phase (degrees). | T |
| RX-08 | The receiver shall be told M; it shall not need any other side information. | I |

---

## 7. Interface requirements

### 7.1 Files

| ID    | Requirement | Verif. |
|-------|-------------|--------|
| IO-01 | The transmitter shall write 16-bit PCM, mono, 48000 Hz WAV files. | T |
| IO-02 | The receiver shall read 16-bit PCM WAV files at 48000 Hz; for multi-channel files it shall use the first channel. | T |
| IO-03 | A WAV with any other sample rate or sample size shall be rejected with an error. | T |
| IO-04 | The C++ harness shall export received symbols as CSV with columns i,q,ideal_i,ideal_q. | T |
| IO-05 | The Python harness shall save a PNG plot (constellation, received symbols, waveform, spectrum). | T |

### 7.2 Command line

```
  qam_demo [--M 4|8|16|32|64] [--text STR] [--snr dB | --clean]
           [--wav out.wav] [--csv rx.csv] [--plot out.png]
           [--rx in.wav] [--sweep]

  --M      constellation order (default 16)
  --text   message to send (default "Hello from QAM!")
  --snr    Es/N0 of the simulated channel in dB (default 30)
  --clean  no channel impairments
  --wav    save the transmitted signal
  --csv    save received symbols (C++ only)
  --plot   save a PNG figure (Python only)
  --rx     decode a WAV file instead of simulating
  --sweep  print a BER / CRC table for all M over a range of SNRs

  Exit codes (C++):  0 = CRC OK    2 = CRC failed    1 = error
```

| ID    | Requirement | Verif. |
|-------|-------------|--------|
| IO-06 | Both harnesses shall accept the options in 7.2 (except the implementation-specific ones marked). | T |
| IO-07 | The C++ harness shall return the exit codes in 7.2 so scripts can test the result. | T |

---

## 8. Performance requirements

Measured with the text "The quick brown fox jumps over the lazy dog",
channel seed 1. Each implementation is tested with its own random
generator (numpy PCG64 in Python, std::mt19937 in C++), so the two see
different noise for the same seed and their results differ slightly.

```
  M     must decode with CRC OK at Es/N0 >=     bit rate
  ---   -----------------------------------     ----------
   4                 12 dB                       4800 bit/s
   8                 15 dB                       7200 bit/s
  16                 18 dB                       9600 bit/s
  32                 24 dB                      12000 bit/s
  64                 24 dB                      14400 bit/s
```

| ID      | Requirement | Verif. |
|---------|-------------|--------|
| PERF-01 | Each M shall decode with CRC OK at the Es/N0 values above and at every higher value up to 30 dB. | T |
| PERF-02 | In clean mode, every M shall decode with zero bit errors. | T |
| PERF-03 | A clean WAV round trip (write, read, decode) shall decode with zero bit errors. | T |
| PERF-04 | The full C++ sweep (5 orders x 9 SNRs) should complete in under 1 second on a laptop. | T |

---

## 9. Code architecture

### 9.1 Module structure (C++)

```
  +-------------------------------------------------------------+
  | qam_demo.cpp   (harness: files, console, simulation)        |
  |                                                             |
  |   main() --+--> run_once() --> channel()                    |
  |            +--> read_wav() / write_wav()                    |
  |            +--> sweep table, CSV export                     |
  +------------------------------+------------------------------+
                                 | uses
                                 v
  +-------------------------------------------------------------+
  | qam_modem.h    (core: no I/O, standard library only)        |
  |   namespace qam                                             |
  |                                                             |
  |   parameters   FS BAUD SPS FC BETA SPAN PREAMBLE_LEN        |
  |                GUARD TX_PEAK                                |
  |   tables       PREAMBLE_IDX  rrc()  preamble()              |
  |   integrity    crc16()                                      |
  |   mapping      class Constellation                          |
  |                  map()  slice()  append_bits()              |
  |   framing      text_to_bits()  bits_to_text()  Frame        |
  |   transmit     modulate()        stages 1-4                 |
  |   receive      demodulate()      stages 5-9 -> RxResult     |
  +-------------------------------------------------------------+
```

### 9.2 Architecture requirements

| ID     | Requirement | Verif. |
|--------|-------------|--------|
| ARC-01 | The modem core shall perform no file, console or network I/O. | I |
| ARC-02 | The C++ core shall depend only on the C++17 standard library. | I |
| ARC-03 | The Python version shall depend only on numpy (matplotlib optional, for plots). | I |
| ARC-04 | Filter taps and the preamble shall be computed once and reused. | I |
| ARC-05 | The core shall be organised top to bottom in signal order, as in 9.1. | I |
| ARC-06 | The C++ harness shall build, from `cpp-poc/`, with: `clang++ -std=c++17 -O2 qam_demo.cpp -o qam_demo`. | T |

---

## 10. Known limitations and Phase 2 requirements

### 10.1 Phase 1 limitations

The signal chain (RRC shaping, matched filter, preamble sync, one-tap
equalizer, Gray-coded QAM) is the standard one used by real modems,
and 2400 baud 16-QAM is in the range of V.32 telephone modems. What is
missing is everything a receiver needs to survive a real link, and
the channel simulator (6.5) does not test for any of it.

```
+------------------------------------+--------------------------------------+
| Limitation                         | Planned improvement                  |
+------------------------------------+--------------------------------------+
| Receiver                           |                                      |
|   Batch processing                 | Streaming, sample-by-sample receiver |
|   Carrier phase measured once, on  | Phase/frequency tracking during the  |
|   the preamble; no frequency       | frame (decision-directed PLL)        |
|   offset correction                |                                      |
|   No symbol timing recovery: one   | Timing recovery loop (e.g. Gardner)  |
|   sample every SPS after sync      |                                      |
|   One-tap equalizer (gain + phase  | Adaptive LMS equalizer, 10-20 taps   |
|   only)                            | (handles echoes)                     |
|   No detection threshold: always   | Threshold on correlation peak;       |
|   takes the strongest peak, even   | several frames per recording         |
|   in noise; one frame per file     |                                      |
|   No automatic gain control        | AGC before the slicer                |
|   Brute-force correlation          | FFT-based correlation                |
|   Receiver must be told M          | Header symbol announcing M           |
| Coding                             |                                      |
|   No error correction: one bit     | FEC (convolutional or Reed-Solomon)  |
|   error loses the whole frame      |                                      |
| Testing                            |                                      |
|   Ideal channel: whole-sample      | Add frequency offset, clock offset,  |
|   delay, flat gain, white noise;   | fractional delay and echoes to the   |
|   phase is always a multiple of    | channel simulator                    |
|   45 degrees                       |                                      |
|   BER from one ~376-bit frame      | Many frames or a long random         |
|   (one error = 0.0027)             | payload per SNR point                |
| Signal                             |                                      |
|   Clean-mode MER tops out at       | Longer SPAN if more headroom is      |
|   ~35 dB (8-symbol RRC truncation) | needed                               |
|   Band 4380-7620 Hz will not pass  | Optional low-carrier mode for        |
|   a phone line (300-3400 Hz)       | voice channels                       |
|   Audio band only; not an RF SDR   | I/Q input from an SDR (RTL-SDR,      |
|                                    | PlutoSDR), with large-offset         |
|                                    | carrier recovery                     |
+------------------------------------+--------------------------------------+
```

### 10.2 Measured robustness

Python receiver, 129-byte frame, Es/N0 30 dB, each impairment applied
on its own (one-off test, not part of the harness). "FAIL" = CRC
failed; the number is MER in dB.

```
  Impairment                 Typical real value         4-QAM     16-QAM    64-QAM
  -------------------------  -------------------------  --------  --------  --------
  none (baseline)            -                          OK 33.4   OK 32.9   OK 33.5
  carrier offset 1 Hz        ~0 via sound card;         OK 33.4   OK 32.9   OK 33.5
  carrier offset 10 Hz       100s of Hz to kHz on RF    FAIL 6.6  FAIL 10.0 FAIL 15.1
  clock offset 20 ppm        10-100 ppm between two     OK 20.0   OK 25.2   OK 27.5
  clock offset 100 ppm       sound cards                FAIL 6.2  FAIL 13.2 FAIL 18.4
  fractional delay 0.5 smp   always present             OK 27.8   OK 28.5   OK 27.0
  echo 0.2 ms, 30%           cables, filters            OK 19.5   OK 19.8   FAIL 18.2
  echo 1 ms, 30%             room acoustics             OK 11.8   FAIL 12.4 FAIL 17.0
  echo 5 ms, 30%             room acoustics             OK  9.9   FAIL 9.9  FAIL 14.1
```

Carrier offsets of about 3 Hz already break 4- and 16-QAM. Clock
drift grows with frame length, so longer messages fail at smaller
offsets than shown.

Expected behaviour on real links today:

```
  Link                                  Expected result
  ------------------------------------  -------------------------------------------
  WAV file written by either harness    Decodes
  Sound-card cable loopback             Probably, for short frames at 4/16-QAM
  Speaker to microphone across a room   Unlikely beyond 4-QAM (echoes)
  Telephone line                        No (signal is outside the voice band)
```

### 10.3 To-do

In priority order:

1. Add frequency offset, clock offset, fractional delay and echoes to
   the channel simulator (CH-*), and include them in `--sweep`, so
   these failures show up in tests before they show up on hardware.
2. Track carrier phase and frequency during the frame (decision-directed
   PLL).
3. Recover symbol timing (e.g. Gardner loop).
4. Replace the one-tap equalizer with an adaptive LMS equalizer
   (10-20 taps).
5. Add forward error correction.
6. Add a detection threshold, multiple frames per recording and AGC.
7. Measure BER over enough bits for a meaningful curve.
8. Announce M in a header symbol so the receiver no longer needs it.

Items 2-4 are also required for Phase 2 (P2-04).

### 10.4 Phase 2 (Arduino) - draft, to be confirmed

| ID    | Requirement | Verif. |
|-------|-------------|--------|
| P2-01 | [P2] The core logic shall be ported to a streaming, fixed-point form suitable for a microcontroller. | T |
| P2-02 | [P2] Filter taps, constellations and the preamble shall be stored as constant tables in flash. | I |
| P2-03 | [P2] A board-generated signal, captured on a computer, shall decode with the Phase 1 receiver (cross-check). | T |
| P2-04 | [P2] The receiver shall recover carrier frequency offset and symbol timing over a real audio link. | T |
| P2-05 | [P2] The target board shall be selected (e.g. Uno R4, Due, ESP32, Teensy 4) based on DAC/ADC and RAM needs. | A |

---

## 11. Verification summary

`qam_demo` below means either `cpp-poc/qam_demo` (built per ARC-06) or
`python3 py-poc/qam_demo.py`.

```
  Test                                   Covers
  -------------------------------------  -------------------------------
  qam_demo --M <m> --clean   (each m)    SYS-03 MAP-* TX-* RX-* PERF-02
  qam_demo --sweep                       PERF-01 PERF-04 CH-*
  C++ --wav -> Python --rx               SYS-04 IO-01 IO-02 PERF-03
  Python --wav -> C++ --rx               SYS-04 IO-01 IO-02 PERF-03
  qam_demo --M 12                        SYS-03 (rejected)
  qam_demo --csv / --plot                IO-04 IO-05
  exit code check in a script            IO-07
```
