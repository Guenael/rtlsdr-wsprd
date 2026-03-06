# WSPR Decoder (`wsprd/`) -- Full Technical Report

## Context

This report provides a comprehensive analysis of the WSPR (Weak Signal Propagation Reporter) decoder implemented in the `wsprd/` directory. WSPR is an amateur radio protocol designed by Joe Taylor (K1JT) to detect extremely weak signals -- up to 28 dB below the noise floor -- using 2-minute transmissions in a 200 Hz bandwidth. The decoder takes 120 seconds of baseband I/Q samples (375 sps) and extracts callsign, grid locator, and transmit power from one or more overlapping WSPR transmissions.

---

## 1. High-Level Architecture

```
+---------------------------------------------------------------------+
|                     WSPR DECODE PIPELINE                            |
|                                                                     |
|  I/Q Samples (45000 @ 375 sps, 120s)                               |
|       |                                                             |
|       v                                                             |
|  +-------------+    +--------------+    +--------------------+      |
|  | FFT Spectrum |--->|  Candidate   |--->|  Coarse Sync       |      |
|  |  Analysis    |    |  Detection   |    |  (FFT grid search) |      |
|  +-------------+    +--------------+    +--------+-----------+      |
|                                                   |                  |
|       +-------------------------------------------+                  |
|       v                                                             |
|  +-------------+    +--------------+    +--------------------+      |
|  | Fine Time   |--->| Fine Freq    |--->| Soft Symbol        |      |
|  | Refinement  |    | Refinement   |    | Extraction + Jitter|      |
|  | (mode 0)    |    | (mode 1)     |    | (mode 2)           |      |
|  +-------------+    +--------------+    +--------+-----------+      |
|                                                   |                  |
|       +-------------------------------------------+                  |
|       v                                                             |
|  +-------------+    +--------------+    +--------------------+      |
|  | De-         |--->| Fano         |--->| Message Unpack     |      |
|  | interleave  |    | Decoder      |    | (bits->call/grid/  |      |
|  |             |    | (K=32, r=1/2)|    |  power)            |      |
|  +-------------+    +--------------+    +--------+-----------+      |
|                                                   |                  |
|       +-------------------------------------------+                  |
|       v                                                             |
|  +-----------------+    +--------------------------+                |
|  | Deduplication   |--->| Signal Subtraction        |               |
|  | & Validation    |    | (remove decoded, repeat)  |               |
|  +-----------------+    +--------------------------+                |
|                                                                     |
|  Output: Array of decoded spots (freq, SNR, dt, drift, message)     |
+---------------------------------------------------------------------+
```

### Multi-Pass Strategy

The decoder runs up to 3 passes over the same I/Q buffer:
- **Pass 1-2**: `maxdrift=4` Hz, `minsync2=0.12` -- finds strong/medium signals
- **Pass 3**: `maxdrift=0`, `minsync2=0.10` -- finds weak signals on cleaned buffer

After each successful decode, the signal is coherently subtracted from the buffer, enabling detection of weaker overlapping signals. This iterative subtraction can decode 10+ simultaneous transmissions.

---

## 2. File-by-File Detailed Analysis

### 2.1 `wsprd.h` -- Decoder Interface

**Constants:**

| Constant | Value | Meaning |
|----------|-------|---------|
| `HASHTAB_SIZE` | 32768 | Callsign hash table buckets |
| `FFT_SIZE` | 512 | FFT length (0.73 Hz/bin at 375 sps) |
| `MAX_CANDIDATES` | 200 | Max frequency candidates per pass |
| `MAX_UNIQUES` | 100 | Max unique decoded messages total |

**Key Data Structures:**

- **`struct cand`** -- Candidate signal: `{freq, snr, shift, drift, sync}`
- **`struct decoder_options`** -- Config: `{freq, rcall, rloc, quickmode, npasses, subtraction}`
- **`struct decoder_results`** -- Decoded spot: `{freq, sync, snr, dt, drift, jitter, message, call, loc, pwr, cycles}`

**Function Prototypes:**

| Function | Purpose |
|----------|---------|
| `sync_and_demodulate()` | 3-mode sync/demod engine (time search, freq search, symbol extraction) |
| `subtract_signal()` | Basic per-symbol signal subtraction |
| `subtract_signal2()` | Enhanced coherent subtraction with LPF smoothing |
| `wspr_decode()` | Main entry point -- full decode pipeline |

---

### 2.2 `wsprd.c` -- Core Decoder (855 lines)

This is the heart of the decoder. It implements the entire signal processing chain.

#### 2.2.1 Signal Parameters

```
Duration:     120 seconds
Sample rate:  375 sps (post CIC decimation from 2.4 Msps)
Samples:      45,000
Symbols:      162 (NSYM)
Info bits:    81 (NBITS)
Samples/sym:  256 (NSPERSYM) -> symbol rate = 375/256 = 1.4648 baud
Tone spacing: DF = 375/256 = 1.4648 Hz
Time/sample:  DT = 1/375 = 2.667 ms
```

WSPR uses **4-FSK** modulation with tones spaced at `DF = 1.4648 Hz`. Each of the 162 symbols selects one of 4 tones (values 0-3). The occupied bandwidth is approximately 6 Hz.

#### 2.2.2 `sync_and_demodulate()` -- The Core Engine

This function operates in three modes, all sharing the same matched-filter architecture:

**Matched Filter Principle:**
For each symbol, four complex correlators run in parallel -- one per tone frequency. Each correlator multiplies the received I/Q signal by a rotating phasor at the tone frequency and integrates over `NSPERSYM=256` samples:

```
z_k[i] = SUM_{j=0}^{255} [ (I[n] + jQ[n]) * exp(-j*2*pi*f_k*j*DT) ]

where f_k = f0 + drift_correction + (k - 1.5) * DF,  k = 0,1,2,3
```

The drift correction models linear frequency drift centered at symbol 81:
```
fp = f0 + (drift/2) * (i - NBITS) / NBITS
```

**Sync Metric Computation:**
The synchronization quality is measured using the `pr3vector` pseudo-random pattern:

```
cmet = (|z_1| + |z_3|) - (|z_0| + |z_2|)     // odd-even tone power difference
ss += (pr3vector[i] == 1) ? +cmet : -cmet      // accumulate with PN sign
sync = ss / total_power                         // normalize to [0, 1]
```

This exploits the fact that the WSPR sync pattern creates a predictable odd/even tone bias. A correctly-aligned signal produces a large positive sync value; noise produces ~0.

**Mode 0 -- Coarse Time Search:**
- Fixes frequency, sweeps time lag from `lagmin` to `lagmax` (step `lagstep`)
- Returns best `shift` (sample offset) and corresponding sync value
- Typical range: +/-128 samples around initial estimate, step 8

**Mode 1 -- Fine Frequency Search:**
- Fixes time lag, sweeps frequency from `ifmin` to `ifmax` (step `fstep=0.1 Hz`)
- Returns best `freq` and updated sync
- Typical range: +/-0.2 Hz around coarse estimate

**Mode 2 -- Soft Symbol Extraction:**
- Fixes both time and frequency
- Extracts soft-decision symbols from tone power differences:
  ```
  fsymb[i] = (pr3[i]==1) ? (|z_3| - |z_1|) : (|z_2| - |z_0|)
  ```
- Normalizes to zero-mean unit-variance, scales to 8-bit range [0, 255]
- `symfac=50` controls scaling factor

**Performance Optimization:** Sin/cos values are cached (`fplast`) and only recomputed when the frequency changes. The Chebyshev recursion `c[j] = c[j-1]*cos(dphi) - s[j-1]*sin(dphi)` avoids per-sample trig calls.

#### 2.2.3 FFT-Based Candidate Search

**Step 1 -- Short-Time FFT:**
- 512-point FFT with Hann window, stepped by 128 samples (quarter-symbol)
- Power spectrum: `ps[bin][block] = |FFT[bin]|^2`
- Frequency resolution: 375/512 = 0.73 Hz/bin (half the tone spacing)

**Step 2 -- Noise Floor Estimation:**
- Average power spectrum across all time blocks
- Smooth with 7-point window
- 30th percentile of smoothed spectrum = noise floor

**Step 3 -- SNR & Peak Detection:**
- `snr[bin] = 10*log10(ps_avg[bin]/noise - 1) - 26.3 dB`
- Local maxima in smoothed spectrum become candidates
- Filtered to +/-110 Hz bandwidth, sorted by SNR descending

#### 2.2.4 Coarse Parameter Estimation via FFT Grid

For each candidate frequency, a 3D grid search over the FFT power array:
- **Frequency**: +/-1 bin around candidate
- **Time**: +/-10 symbols (+/-21 FFT blocks)
- **Drift**: +/-`maxdrift` bins

The sync metric is computed using the same `pr3vector` pattern but directly on FFT bin amplitudes (no demodulation). This provides rough `{shift, drift, freq}` estimates.

#### 2.2.5 Iterative Refinement and Decoding

For each candidate passing the sync threshold (`minsync1=0.10`):

1. **Mode 0**: Refine time lag (+/-128 samples, step 8)
2. **Mode 1**: Refine frequency (+/-0.2 Hz, step 0.1 Hz)
3. **Jittered decode loop** (up to ~42 attempts):
   - Apply time jitter: 0, +/-3, +/-6, ... +/-128 samples
   - **Mode 2**: Extract soft symbols at jittered time
   - Check symbol RMS > `minrms=52` and sync > `minsync2=0.12`
   - **De-interleave** symbols
   - **Fano decode** with `delta=60`, `maxcycles=10000`
   - If successful: unpack message, validate, deduplicate

#### 2.2.6 Signal Subtraction

**`subtract_signal()`** -- Basic method:
- For each symbol, correlate received signal with reference tone to extract complex amplitude
- Regenerate ideal signal with that amplitude and subtract sample-by-sample
- Simple but leaves residual noise due to per-symbol amplitude discontinuities

**`subtract_signal2()`** -- Enhanced method (used in practice):
- Generates continuous-phase reference signal for all 162 symbols
- Computes complex baseband amplitude: `c(t) = s(t) * conj(r(t))`
- Low-pass filters with 360-tap half-sine window (smooths amplitude estimate)
- Edge correction normalizes partial window responses at signal boundaries
- Subtracts: `s'(t) = s(t) - LPF[c(t)] * r(t)`
- Produces cleaner subtraction with fewer artifacts

#### 2.2.7 Key Thresholds

| Parameter | Value | Purpose |
|-----------|-------|---------|
| `minsync1` | 0.10 | FFT-pass sync threshold |
| `minsync2` | 0.12 | Fano-attempt sync threshold |
| `minrms` | 52.0 | Minimum symbol amplitude |
| `maxdrift` | 4 | Max frequency drift (Hz) |
| `symfac` | 50 | Soft symbol scaling |
| `delta` | 60 | Fano threshold step |
| `maxcycles` | 10000 | Fano iteration limit |
| `fmin/fmax` | +/-110 Hz | Candidate frequency range |

---

### 2.3 `fano.c` / `fano.h` -- Fano Sequential Decoder (238 lines)

#### 2.3.1 Convolutional Code Parameters

| Parameter | Value |
|-----------|-------|
| Rate | 1/2 (1 input bit -> 2 output symbols) |
| Constraint length K | 32 (32-bit shift register) |
| Polynomials | Layland-Lushbaugh: POLY1=`0xf2d05351`, POLY2=`0xe4613c47` |
| Code type | Nonsystematic, non-quick-look-in |

Two alternative polynomial sets are available via preprocessor (NASA Standard, Massey-Johannesson) but the active code uses Layland-Lushbaugh (LL).

#### 2.3.2 The ENCODE Macro

Computes 2-bit output from encoder state:
```
sym[bit1] = parity(encstate & POLY1)    // XOR of selected state bits
sym[bit0] = parity(encstate & POLY2)
```
Parity computed via cascaded XOR-shifts and `Partab[]` lookup (256-entry table in `tab.c`).

#### 2.3.3 Node Structure

```c
struct node {
    unsigned long encstate;   // 32-bit encoder shift register state
    long gamma;               // Cumulative path metric (cost to reach this node)
    int metrics[4];           // Branch metrics for all 4 possible symbol pairs
    int tm[2];                // Sorted: tm[0]=better branch, tm[1]=worse branch
    int i;                    // Current branch index (0=better, 1=worse)
};
```

#### 2.3.4 Fano Algorithm

The Fano algorithm is a **sequential tree-search decoder** that explores the code trellis one branch at a time, using a dynamic threshold to prune unlikely paths:

**Initialization:**
- Allocate `nbits+1` nodes (one per decoded bit)
- **Precompute all branch metrics** from received symbols -- this is the only time raw symbols are accessed:
  ```
  metrics[0] = mettab[0][sym0] + mettab[0][sym1]   // both bits = 0
  metrics[1] = mettab[0][sym0] + mettab[1][sym1]   // bit0=0, bit1=1
  metrics[2] = mettab[1][sym0] + mettab[0][sym1]   // bit0=1, bit1=0
  metrics[3] = mettab[1][sym0] + mettab[1][sym1]   // both bits = 1
  ```
- Start at root with `gamma=0`, `threshold=0`

**Main Loop (one cycle per iteration, bounded by `maxcycles`):**

```
FORWARD: Try current best branch
  ngamma = gamma + tm[i]
  IF ngamma >= threshold:
    Tighten threshold (while ngamma >= threshold + delta: threshold += delta)
    Move forward to next node
    Compute new branch metrics, sort into tm[0] (better) / tm[1] (worse)
    IF reached last node -> SUCCESS

BACKWARD: Current branch failed threshold
  Look back through previous nodes:
    IF parent's gamma < threshold:
      Lower threshold (threshold -= delta)
      Reset to better branch
    ELSE IF untried branch exists:
      Try alternate branch
    ELSE:
      Continue backtracking
```

**Tail Constraint:** After bit `nbits-31`, all remaining bits are forced to 0 (the K-1 = 31 zero-tail flushes the encoder state). Only the 0-branch is explored in the tail region.

**Output:** Decoded bits extracted from encoder states at 8-bit intervals. Returns 0 on success, -1 on timeout.

**Why Fano instead of Viterbi?** For K=32, Viterbi would require 2^31 = 2 billion states -- computationally infeasible. Fano explores only the most likely path, with occasional backtracking, making it practical for very long constraint lengths.

#### 2.3.5 `encode()` Function

Reference convolutional encoder used for self-test and signal subtraction:
- Iterates through input bytes MSB-first
- Shifts each bit into encoder state
- Applies ENCODE macro to produce 2 symbols per bit
- Outputs `2 * 8 * nbytes` symbols

---

### 2.4 `wsprd_utils.c` / `wsprd_utils.h` -- Message Unpacking (314 lines)

These files handle the **decoder side**: converting 81 decoded bits back into human-readable callsign, grid locator, and power.

#### 2.4.1 `unpack50()` -- Bit Extraction

Extracts 50 information bits from 11-byte Fano output into two integers:
- **n1** (28 bits): callsign encoding
- **n2** (22 bits): grid + power + message type

```
Byte layout: [dat0][dat1][dat2][dat3][dat4][dat5][dat6][dat7-10=tail]
              |-- n1 (28 bits) --|---- n2 (22 bits) ---|
```

#### 2.4.2 `unpackcall()` -- Callsign Decoding

Decodes 28-bit integer to 6-character callsign using **mixed-radix** decomposition:

```
Position:  [0]    [1]    [2]    [3]    [4]    [5]
Base:       36     36     10     27     27     27
Charset:   0-9,   0-9,   0-9,   A-Z,   A-Z,   A-Z,
           A-Z    A-Z           space  space  space

Max value: 36 x 36 x 10 x 27 x 27 x 27 = 262,177,560 (fits in 28 bits)
```

This encoding exploits the structure of amateur callsigns: position 2 must be a digit, positions 3-5 are letters or space.

#### 2.4.3 `unpackgrid()` -- Grid Locator Decoding

Decodes 18 bits (from n2 >> 7) to 4-character Maidenhead grid:

```
ngrid < 32400
dlat  = (ngrid % 180) - 90          -> latitude
dlong = (ngrid / 180) * 2 - 180     -> longitude

Grid[0] = longitude field (A-R)
Grid[1] = latitude field (A-R)
Grid[2] = longitude square (0-9)
Grid[3] = latitude square (0-9)
```

#### 2.4.4 `unpackpfx()` -- Prefix/Suffix Handling

For non-standard callsigns (e.g., `PJ4/K1ABC`, `W1AW/7`):
- **nprefix < 60000**: Prefix mode -- decode 3 base-37 digits -> "PFX/CALL"
- **nprefix >= 60000**: Suffix mode -- decode 1-2 characters -> "CALL/SFX"

#### 2.4.5 `unpk_()` -- Master Unpacking

Orchestrates the full unpack pipeline. Determines message type from `ntype = (n2 & 127) - 64`:

| Type | Condition | Format | Grid Precision |
|------|-----------|--------|----------------|
| Type 1 (Standard) | `ntype` in {0,3,7,10,13,...,60} | `CALL GRID4 POWER` | 4-char (~2.4 km) |
| Type 2 (Extended) | `ntype > 0`, not standard power | `PFX/CALL POWER` | None (traded for prefix) |
| Type 3 (Hashed) | `ntype < 0` | `<HASH> GRID6 POWER` | 6-char (~200 m) |

Type 3 messages sacrifice callsign information (replaced by hash) to transmit a 6-character grid for higher location precision. The hash table allows recovery of known callsigns.

#### 2.4.6 `deinterleave()` -- Reverse Permutation

Reverses the channel interleaving using a bit-reversal permutation:
```c
j = ((i * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> 32
```
This elegant constant-multiplication trick computes the bit-reversal of an 8-bit index without loops.

---

### 2.5 `wsprsim_utils.c` / `wsprsim_utils.h` -- Message Packing & Encoding (317 lines)

These files handle the **encoder side**: converting callsign/grid/power into 162 channel symbols. Used for self-test and signal subtraction (re-encoding decoded messages).

#### 2.5.1 `pack_call()` -- Callsign Encoding

Inverse of `unpackcall()`: 6-char callsign -> 28-bit integer via mixed-radix:
```
n = c[0]*36*10*27*27*27 + c[1]*10*27*27*27 + c[2]*27*27*27 + c[3]*27*27 + c[4]*27 + c[5]
```

#### 2.5.2 `pack_grid4_power()` -- Grid+Power Encoding

Encodes 4-char grid + power to 22-bit value:
```
m = (179 - 10*grid[0] - grid[2]) * 180 + 10*grid[1] + grid[3]
m = m * 128 + power + 64
```

#### 2.5.3 `get_wspr_channel_symbols()` -- Complete Encoding Pipeline

This is the full encoder, producing 162 channel symbols from a text message:

```
Step 1: Parse message -> determine type (1/2/3)
Step 2: Pack callsign -> n (28 bits)
Step 3: Pack grid+power -> m (22 bits)
Step 4: Pack n,m into 11 bytes (50 info bits + 31-bit zero tail)
Step 5: Convolutional encode -> 162 raw symbols (rate 1/2, K=32)
Step 6: Interleave (bit-reversal permutation)
Step 7: Merge with sync pattern:
        channel_symbols[i] = 2 * data[i] + pr3vector[i]
        -> produces values 0, 1, 2, 3 (4-FSK tones)
```

**Power Level Quantization:** WSPR power must be from the set {0, 3, 7, 10, 13, 17, 20, 23, 27, 30, 33, 37, 40, 43, 47, 50, 53, 57, 60} dBm. An adjustment table snaps arbitrary values to the nearest valid level:
```
nu[10] = {0, -1, 1, 0, -1, 2, 1, 0, -1, 1}
power += nu[power % 10]
```

#### 2.5.4 `interleave()` -- Channel Interleaving

Forward permutation (inverse of `deinterleave()`). Spreads adjacent encoded bits across the 162-symbol frame, protecting against burst errors and short fades on HF channels.

---

### 2.6 `nhash.c` / `nhash.h` -- Jenkins Hash (451 lines)

Implements Bob Jenkins' **lookup3** hash algorithm for callsign hashing.

**Purpose in WSPR:** Type 3 messages encode a 15-bit hash of the callsign instead of the full callsign. The receiver maintains a hash table of previously-heard callsigns to resolve hashes back to callsigns.

**Algorithm:**
1. Initialize 3 accumulators: `a = b = c = 0xdeadbeef + length + initval`
2. Process 12 bytes per round with `mix(a,b,c)` -- 6 reversible operations
3. Handle tail (0-11 bytes) with masking
4. Apply `final(a,b,c)` -- 7-step avalanche mixing
5. Mask to 15 bits: `c = c & 0x7FFF`

Three code paths for different memory alignments (32-bit, 16-bit, byte-by-byte) optimize performance. The hash is non-cryptographic but provides excellent distribution for table lookups.

---

### 2.7 `tab.c` -- Parity Lookup Table (42 lines)

256-entry `Partab[]` array: `Partab[i]` = parity (popcount mod 2) of byte `i`.

Used by the ENCODE macro in `fano.h` to compute convolutional encoder output bits efficiently. The cascaded XOR-shift `(_tmp ^ (_tmp >> 16) ^ (_tmp >> 8)) & 0xFF` reduces a 32-bit parity computation to a single table lookup.

---

### 2.8 `metric_tables.h` -- Soft Decision Metrics (139 lines)

Pre-computed `float metric_tables[5][256]` -- log-likelihood ratio tables for 5 SNR levels:

| Index | Es/No | Use Case |
|-------|-------|----------|
| 0 | 0 dB | Very weak signals |
| 1 | 3 dB | Weak signals |
| 2 | 6 dB | **Typical (default)** |
| 3 | 9 dB | Strong signals |
| 4 | 12 dB | Very strong signals |

**Usage in decoder:**
```c
bias = 0.45
mettab[0][i] = round(10 * (metric_tables[2][i] - bias))     // P(bit=0 | symbol=i)
mettab[1][i] = round(10 * (metric_tables[2][255-i] - bias))  // P(bit=1 | symbol=i)
```

The Es/No=6dB table (index 2) is used with a 0.45 bias offset. The reversal `255-i` for `mettab[1]` exploits the symmetry of the soft-decision metric.

---

## 3. Algorithm Summary

### 3.1 WSPR Modulation

```
Message: "K1ABC EN50 33"
    | pack (50 bits)
[28-bit callsign | 22-bit grid+power]
    | zero-pad (31-bit tail)
81 bits total
    | convolutional encode (K=32, r=1/2)
162 binary symbols
    | interleave (bit-reversal permutation)
162 interleaved symbols
    | merge with pr3 sync pattern
162 quaternary symbols (0-3)
    | 4-FSK modulate (1.4648 Hz spacing)
~6 Hz bandwidth, 110.6 second transmission
```

### 3.2 WSPR Demodulation

```
45000 I/Q samples (120s @ 375 sps)
    | Hann-windowed 512-pt FFT (stepped by 128 samples)
Power spectrum [411 bins x ~320 time blocks]
    | noise floor estimation (30th percentile)
    | local maxima detection
~20-50 candidates (freq, SNR)
    | coarse sync (3D grid: freq +/-1 bin, time +/-10 sym, drift +/-4 Hz)
Rough {shift, drift, freq} per candidate
    | fine time refinement (sync_and_demodulate mode 0)
    | fine freq refinement (sync_and_demodulate mode 1)
Precise {shift, freq} per candidate
    | jittered symbol extraction (mode 2, up to ~42 offsets)
    | deinterleave
    | Fano sequential decode (K=32, delta=60, maxcycles=10000)
81 decoded bits
    | unpack (mixed-radix callsign + grid/power decomposition)
    | deduplicate and validate
Decoded spot: {call, grid, power, freq, SNR, dt, drift}
    | subtract_signal2 (coherent LPF subtraction)
Cleaned buffer -> repeat for next pass
```

### 3.3 Why It Works at -28 dB SNR

The WSPR protocol achieves extreme sensitivity through:

1. **Long integration** -- 120-second symbols provide 50 dB processing gain over 1-second bandwidth
2. **Strong FEC** -- K=32 convolutional code with Fano soft-decision decoding approaches Shannon limit
3. **Narrow bandwidth** -- 6 Hz total, 1.46 Hz tone spacing minimizes noise
4. **Coherent sync** -- pr3 pseudo-random pattern enables robust synchronization even below noise
5. **Iterative subtraction** -- removes strong signals to reveal weaker ones beneath
6. **Soft decisions** -- 8-bit symbol quantization preserves analog confidence for the decoder
7. **Minimal payload** -- only 50 bits of information (callsign + grid + power) maximizes coding gain

---

## 4. File Summary Table

| File | Lines | Role |
|------|-------|------|
| `wsprd.h` | 111 | Public API: structs, constants, function prototypes |
| `wsprd.c` | 855 | Core decoder: FFT search, sync, demod, subtraction, orchestration |
| `fano.h` | 45 | Fano decoder API + ENCODE macro + polynomial definitions |
| `fano.c` | 238 | Fano sequential decoder + convolutional encoder |
| `wsprd_utils.h` | 43 | Unpacking function prototypes |
| `wsprd_utils.c` | 314 | Message unpacking: bits -> callsign/grid/power |
| `wsprsim_utils.h` | 10 | Packing/encoding function prototypes |
| `wsprsim_utils.c` | 317 | Message packing + full encoding pipeline (used for self-test & subtraction) |
| `nhash.h` | 3 | Hash function prototype |
| `nhash.c` | 451 | Bob Jenkins lookup3 hash (15-bit callsign hashing) |
| `metric_tables.h` | 139 | Pre-computed soft-decision LLR tables (5 SNR levels x 256 entries) |
| `tab.c` | 42 | 256-byte parity lookup table for ENCODE macro |
