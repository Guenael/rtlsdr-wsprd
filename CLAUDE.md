# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

rtlsdr-wsprd is a C daemon that receives and decodes WSPR (Weak Signal Propagation Reporter) amateur radio signals using an RTL-SDR USB dongle, then reports decoded spots to wsprnet.org via HTTPS. Designed to run headless on Raspberry Pi.

## Build Commands

```bash
make                # Build rtlsdr_wsprd (uses clang, -O3, -std=gnu17)
sudo make install   # Install to /usr/local/bin/rtlsdr_wsprd
make clean          # Remove build artifacts
```

Dependencies: `librtlsdr`, `libusb-1.0`, `libfftw3f`, `libcurl`, `libpthread`

On Ubuntu/Debian: `apt-get install librtlsdr-dev libfftw3-dev libusb-1.0-0-dev libcurl4-gnutls-dev`

## Testing

No test framework — tests use the built-in self-test mode:

```bash
# Reference IQ file decode test (requires built binary)
./rtlsdr_wsprd -f 2m -c A1XYZ -l AB12CD -r ./signals/refSignalSnr0dB.iq

# Synthetic signal self-test
./rtlsdr_wsprd -f 2m -c A1XYZ -l AB12CD -t
```

Both return exit code 0 on success, 1 on failure.

## CI

GitHub Actions (`.github/workflows/ci.yml`): CodeQL, cppcheck, cpplint, codespell, then both unit tests above.

## Architecture

### Signal Processing Pipeline

`rtlsdr_wsprd.c` → `wsprd/wsprd.c`

1. **RTL-SDR capture** at 2.4 Msps (8-bit unsigned IQ)
2. **fs/4 mixer** — shifts tuned frequency by 600 kHz using sign-flip trick (no trig)
3. **CIC decimation** (N=2, order=2) — 2.4 Msps down to 375 sps (factor 6400)
4. **FIR compensation** (32-tap symmetric) — corrects CIC frequency droop
5. **WSPR decode** — FFT candidate search → coarse/fine sync → Fano convolutional decoder (K=32) → message unpack
6. **HTTP POST** decoded spots to wsprnet.org

### Threading Model (3 threads)

- **Main thread**: Initializes RTL-SDR, runs 2-minute timing loop aligned to UTC even minutes, flips double buffer, signals decoder
- **RX thread** (`rtlsdr_rx`): Runs `rtlsdr_read_async()`, callback does mixing + CIC decimation into current buffer
- **Decoder thread**: Waits on condition variable, processes previous buffer, calls `wspr_decode()`, posts results

Synchronization via `safe_cond_signal`/`safe_cond_wait` macros (mutex + condvar). Decoder runs at SCHED_RR priority 90.

### Double Buffering

`rx_state.iSamples[2][45000]` / `qSamples[2][45000]` — 120s × 375 sps. `bufferIndex` alternates 0/1 every 2 minutes so RX and decode run in parallel on separate buffers.

### Key Source Files

- `rtlsdr_wsprd.c` — Main app: SDR I/O, threading, CLI args, HTTP spot reporting, CIC/FIR DSP
- `wsprd/wsprd.c` — WSPR decoder: FFT spectrum analysis, sync, soft-symbol extraction, Fano decode, iterative subtraction
- `wsprd/wsprd_utils.c` — Message unpacking (callsign/grid/power from decoded bits)
- `wsprd/wsprsim_utils.c` — Message packing / channel symbol generation (encoder side, used by self-test)
- `wsprd/fano.c` — Fano sequential convolutional decoder (KA9Q)
- `wsprd/nhash.c` — Jenkins hash for callsign hash table

## Code Style

Google-based clang-format (`.clang-format`): 4-space indent, no tabs, `clang` compiler.
