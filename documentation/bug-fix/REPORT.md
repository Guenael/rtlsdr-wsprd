# Rapport: rtlsdr-wsprd Code Review & Fixes

Implementation of all four phases from `PLAN.md`. All changes verified against the self-test, reference IQ decode, and `--version` exit code.

---

## Phase 1: Bug Fixes

### 1A. Self-test logic inversion — `rtlsdr_wsprd.c`

**Bug:** `strcmp` returns 0 on match, so `&&` between three `strcmp` calls made the condition true only when all three were non-zero (all mismatches). Due to short-circuit evaluation, a single matching field (returning 0) would short-circuit the whole expression to false, incorrectly reporting success.

**Fix:** Changed `&&` to `||` so the test correctly fails when at least one field doesn't match.

**Status:** Done, verified by self-test.

### 1B. `--version` exit code & missing break — `rtlsdr_wsprd.c`

**Bug:** `exit(EXIT_FAILURE)` after printing version. Also, `case 0:` (long options) was missing a `break`, causing fallthrough into `case 'f':` (frequency parsing).

**Fix:** Changed to `exit(EXIT_SUCCESS)`, added `break` after the inner switch.

**Status:** Done, `--version` now exits 0.

### 1C. Bitwise `&` → logical `&&` — `wsprd.c`, `wsprd_utils.c`

**Bug:** 7 instances of `if ((x >= 0) & (x < N))` using bitwise AND instead of logical AND. While functionally equivalent for boolean operands, this is undefined-behavior-adjacent and defeats short-circuit evaluation.

**Fix:** Changed all 7 instances to `&&`.

**Files changed:**
- `wsprd/wsprd.c` — 2 instances (lines 281, 293)
- `wsprd/wsprd_utils.c` — 5 instances (lines 165, 167, 183, 188, 193)

**Status:** Done.

---

## Phase 2: Safety Hardening

### 2A. Hash table index validation — `wsprd/wsprd.c`

**Bug:** `sscanf` reads `nh` from `hashtable.txt` and uses it directly as an array index without bounds checking. A malformed file could cause out-of-bounds writes.

**Fix:** Added `if (nh >= 0 && nh < HASHTAB_SIZE)` guard before the `strcpy` calls.

**Status:** Done.

### 2B. fopen NULL check — `wsprd/wsprd.c`

**Bug:** `fopen("hashtable.txt", "w")` result was not checked before `fprintf`. If the file can't be opened (e.g., read-only filesystem), this is a NULL pointer dereference.

**Fix:** Wrapped the write loop in `if (fhash) { ... fclose(fhash); }`.

**Status:** Done.

### 2C. `volatile` on `exit_flag` — `rtlsdr_wsprd.c`

**Bug:** `exit_flag` is set in a signal handler and read in multiple threads, but was not declared `volatile`. The compiler could optimize away the reads.

**Fix:** Changed `bool exit_flag` to `volatile bool exit_flag`.

**Status:** Done.

### 2D. malloc NULL checks & memory leak — `wsprd/wsprsim_utils.c`

**Bug:** 5 `malloc` calls (lines 281–285) had no NULL checks. Additionally, `call`, `loc`, `pwr` (lines 283–285) were never freed — only `check_call_loc_pow` and `check_callsign` were freed.

**Fix:** Added NULL check after all 5 mallocs with early return on failure. Added the 3 missing `free()` calls.

**Status:** Done.

### 2E. URL parameter escaping — `rtlsdr_wsprd.c`

**Bug:** User-supplied `rcall` and `rloc` fields were embedded directly in URLs without escaping. Special characters could break the HTTP request or enable injection.

**Fix:** Used `curl_easy_escape()` for both fields before embedding in URLs. Also refactored `postSpots()` to reuse a single curl handle and increased URL buffer from 256 to 512 bytes.

**Status:** Done.

### 2F. Unsafe `strcpy`/`sprintf` → bounded alternatives

**Bug:** Multiple `strcpy` and `sprintf` calls could overflow their destination buffers:
- `wsprd_utils.c:159` — `strcpy(tmpcall, call)` into `char tmpcall[7]` when call can be 12 chars
- `wsprd_utils.c:325` — `sprintf(callsign, "<%s>", ...)` could overflow 13-byte buffer
- `wsprd_utils.c:274–275,302` — `strcpy` into hashtab/loctab without bounds
- `wsprd.c:799,813–816` — `strcpy` into allcalls/decodes arrays
- Various other `strcpy`/`sprintf` throughout `wsprd_utils.c`

**Fix:** Replaced all with `snprintf` bounded calls. Enlarged `tmpcall` to 13 bytes. Replaced `strncat` chains in `unpk_()` with single `snprintf` calls.

**Status:** Done.

---

## Phase 3: Build System & CI

### 3A. Makefile modernization

**Changes:**
- `CC = clang` → `CC ?= clang` (allows override via `make CC=gcc`)
- Added `-Wall -Wextra -Wno-unused-parameter` to CFLAGS
- Added `$(LDFLAGS)` to the link line for distro/packaging compatibility
- Added `install` to `.PHONY`
- Added aarch64 detection block (uses same `-DRPI23` as armv7)
- Added `-MMD` for automatic header dependency tracking with `-include $(DEPS)`
- Clean target now removes `*.d` files

**Status:** Done.

### 3B. CI workflow fixes

**Changes:**
- `container.yml`: `crazy-max/ghaction-docker-meta@v1` → `docker/metadata-action@v5`
- `ci.yml`: Fixed typo "best best-practice" → "best-practice"
- `ci.yml`: Added comment explaining `|| true` on cpplint

**Status:** Done.

### 3C. Dockerfile & .dockerignore

**Changes:**
- Created `.dockerignore` excluding `.git`, `*.o`, `*.d`, `*.iq`, `art/`, `.github`, `CLAUDE.md`, `PLAN.md`
- Pinned librtlsdr download to release tag `v2.0.2` instead of `master.zip`

**Status:** Done.

---

## Phase 4: Code Quality

### 4A. Replace bubble sort with `qsort`

**Change:** Replaced two O(n²) bubble sorts with `qsort`:
1. Candidate sorting by SNR (descending) in the FFT peak search
2. Decoder results sorting by SNR (descending) before returning

Added two static comparison functions: `cand_snr_desc` and `results_snr_desc`.

**Status:** Done.

### 4B. Extract magic numbers into named constants

**Added to `wsprd/wsprd.h`:**

| Constant | Value | Replaces |
|----------|-------|----------|
| `HASHTAB_SIZE` | 32768 | `32768` |
| `HASHTAB_ENTRY_LEN` | 13 | `13` (hash entry width) |
| `LOCTAB_ENTRY_LEN` | 5 | `5` (loc entry width) |
| `FFT_SIZE` | 512 | `512` (FFT bins) |
| `MAX_CANDIDATES` | 200 | `200` (candidate array size) |
| `MAX_UNIQUES` | 100 | `100` (unique decode limit) |

Updated all references in `wsprd.c`, `wsprd_utils.c`, and `rtlsdr_wsprd.c`.

**Status:** Done.

### 4C. `const` qualifiers on read-only pointer params

**Changes in `rtlsdr_wsprd.h` / `rtlsdr_wsprd.c`:**
- `readRawIQfile(..., const char *filename)`
- `writeRawIQfile(..., const char *filename)`
- `readC2file(..., const char *filename)`
- `decodeRecordedFile(const char *filename)`

**Changes in `wsprd/wsprd.h` / `wsprd/wsprd.c`:**
- `subtract_signal(..., const unsigned char *channel_symbols)`
- `subtract_signal2(..., const unsigned char *channel_symbols)`

**Status:** Done.

### 4D. Clean up stale comment markers

**Removed:**
- `rtlsdr_wsprd.c:174` — `// CHECK -127 alt. possible issue ?`
- `rtlsdr_wsprd.c:190` — `// UPDATE: i+=2 & fix below`
- `rtlsdr_wsprd.c:192` — `// EVAL: option to move sigIn in float here`
- `wsprd/wsprd_utils.c:38` — `// EVAL -- Replace strcpy & strncpy ...` (plus two commented-out pragmas)

**Status:** Done.

---

## Bonus fixes (discovered during implementation)

- **`static const` ordering** — `const static float zCoef[]` → `static const float zCoef[]` (fixes `-Wold-style-declaration`)
- **Uninitialized variable** — `float sync` in coarse search loop was uninitialized; added `= 0.0`
- **Filename buffer too small** — `saveSample()` had `char filename[32]` for a format that can produce up to 40 chars; enlarged to 64

---

## Verification

```
$ make clean && make CC=gcc          # builds OK
$ ./rtlsdr_wsprd -f 2m -c A1XYZ -l AB12CD -t
Spot(0)  22.80   0.01 144.490550  0    K1JT   FN20 20
Self-test SUCCESS!                   # exit 0 ✓

$ ./rtlsdr_wsprd -f 2m -c A1XYZ -l AB12CD -r ./signals/refSignalSnr0dB.iq
Spot :  -0.07   0.01 144.490550  0    K1JT   FN20 20
                                     # exit 0 ✓

$ ./rtlsdr_wsprd --version
rtlsdr_wsprd v0.5.6                 # exit 0 ✓
```

## Files modified

| File | Phases |
|------|--------|
| `rtlsdr_wsprd.c` | 1A, 1B, 2C, 2E, 4B, 4C, 4D |
| `rtlsdr_wsprd.h` | 4C |
| `wsprd/wsprd.c` | 1C, 2A, 2B, 2F, 4A, 4B, 4C |
| `wsprd/wsprd.h` | 4B, 4C |
| `wsprd/wsprd_utils.c` | 1C, 2F, 4B, 4D |
| `wsprd/wsprsim_utils.c` | 2D |
| `Makefile` | 3A |
| `.github/workflows/ci.yml` | 3B |
| `.github/workflows/container.yml` | 3B |
| `Dockerfile` | 3C |

## Files created

| File | Phase |
|------|-------|
| `.dockerignore` | 3C |
| `RAPPORT.md` | — |
