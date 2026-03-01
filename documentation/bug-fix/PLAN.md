# Plan: rtlsdr-wsprd Code Review & Fixes

## Context

The project is a C daemon that decodes WSPR amateur radio signals via RTL-SDR and reports to wsprnet.org. It's described as "old and needing love." A thorough review found real bugs, safety issues, and build system rot. This plan addresses them in phases ordered by severity, each independently committable and testable.

---

## Phase 1: Bug Fixes (high-value, safe, immediate)

### 1A. Self-test logic inversion — `rtlsdr_wsprd.c:773-779`

**Bug:** Uses `&&` but needs `||`. Currently passes if ANY ONE field matches instead of requiring ALL THREE to match (because `strcmp` returns 0 on match, and `0 && x` short-circuits to false).

```c
// BEFORE (buggy)
if (strcmp(dec_results[0].call, "K1JT") &&
    strcmp(dec_results[0].loc,  "FN20") &&
    strcmp(dec_results[0].pwr,  "20")) {
    return 0;  // fail
} else {
    return 1;  // pass
}

// AFTER (correct)
if (strcmp(dec_results[0].call, "K1JT") ||
    strcmp(dec_results[0].loc,  "FN20") ||
    strcmp(dec_results[0].pwr,  "20")) {
    return 0;  // fail: at least one mismatch
} else {
    return 1;  // pass: all match
}
```

### 1B. `--version` exit code — `rtlsdr_wsprd.c:849`

Change `exit(EXIT_FAILURE)` to `exit(EXIT_SUCCESS)`.

### 1C. Bitwise `&` → logical `&&` in conditions

| File | Lines |
|------|-------|
| `wsprd/wsprd.c` | 281, 293 |
| `wsprd/wsprd_utils.c` | 165, 167, 183, 188, 193 |

All are `if ((x >= 0) & (x < N))` patterns — change `&` to `&&`.

---

## Phase 2: Safety Hardening

### 2A. Hash table index validation — `wsprd/wsprd.c:474-476`

`sscanf` reads `nh` from `hashtable.txt` and uses it as array index without bounds check. Add `if (nh >= 0 && nh < 32768)` guard before the `strcpy` calls.

### 2B. fopen NULL check — `wsprd/wsprd.c:848`

`fhash = fopen("hashtable.txt", "w")` is not checked before `fprintf(fhash, ...)`. Wrap the write loop in `if (fhash) { ... fclose(fhash); }`.

### 2C. `volatile` on `exit_flag` — `rtlsdr_wsprd.c:77`

Changed in signal handler (line 249), read in threads. Change `bool exit_flag` to `volatile bool exit_flag`.

### 2D. malloc NULL checks & memory leak — `wsprd/wsprsim_utils.c:281-305`

Add NULL checks after 5 malloc calls (lines 281-285). Also fix memory leak: `call`, `loc`, `pwr` (lines 283-285) are never freed — only `check_call_loc_pow` and `check_callsign` are freed at lines 303-304. Add the 3 missing `free()` calls.

### 2E. URL parameter escaping — `rtlsdr_wsprd.c:366-435`

In `postSpots()`, use `curl_easy_escape()` for user-supplied `rcall` and `rloc` fields before embedding in URLs. The curl handle is already available. Free the escaped strings with `curl_free()`.

### 2F. Unsafe `strcpy`/`sprintf` → bounded alternatives

Key locations (most dangerous first):

| File | Line | Issue |
|------|------|-------|
| `wsprd/wsprd_utils.c:159` | `strcpy(tmpcall, call)` into `char tmpcall[7]` — call can be 12 chars |
| `wsprd/wsprd_utils.c:325` | `sprintf(callsign, "<%s>", hashtab+ihash*13)` — can overflow 13-byte callsign |
| `wsprd/wsprd_utils.c:274-275,302` | `strcpy` into hashtab/loctab without index validation |
| `wsprd/wsprd.c:475-476` | `strcpy` into hashtab/loctab (covered by 2A bounds check) |
| `wsprd/wsprd.c:799,813-816` | `strcpy` into allcalls/decodes arrays |
| `wsprd/wsprd_utils.c:84,109,148,176,185,190,196` | Various `strcpy`/`sprintf` |

Replace with `snprintf(dst, sizeof(dst), "%s", src)` or equivalent bounded calls. For pointer destinations where sizeof isn't available, pass explicit size constants.

---

## Phase 3: Build System & CI

### 3A. Makefile — `Makefile`

- `CC = clang` → `CC ?= clang` (allow override)
- Add `-Wall -Wextra -Wno-unused-parameter` to CFLAGS
- Link line: add `$(LDFLAGS)` → `$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)`
- Add `install` to `.PHONY`
- Add aarch64 detection block (use same `-DRPI23` as armv7)
- Add `-MMD` for header dependency tracking, `-include $(DEPS)`, clean `*.d`

### 3B. CI workflows

- `container.yml:21`: `crazy-max/ghaction-docker-meta@v1` → `docker/metadata-action@v5`
- `ci.yml:47`: Fix typo "best best-practice" → "best-practice"
- `ci.yml:53`: Add comment explaining `|| true` on cpplint

### 3C. Dockerfile

- Create `.dockerignore` (exclude `.git`, `*.o`, `*.iq`, `art/`, `.github`)
- Pin librtlsdr download to a release tag instead of `master.zip`

---

## Phase 4: Code Quality

### 4A. Replace bubble sort with `qsort` — `wsprd/wsprd.c:617-627`

Add a `cand_snr_desc` comparison function and use `qsort(candidates, npk, sizeof(struct cand), cand_snr_desc)`.

### 4B. Extract magic numbers into named constants in `wsprd/wsprd.h`

`32768` → `HASHTAB_SIZE`, `13` → `HASHTAB_ENTRY_LEN`, `5` → `LOCTAB_ENTRY_LEN`, `512` → `FFT_SIZE`, `200` → `MAX_CANDIDATES`, `100` → `MAX_UNIQUES`

### 4C. Add `const` qualifiers to read-only pointer params in headers

- `rtlsdr_wsprd.h`: `atofs(const char *s)`, `readRawIQfile(..., const char *filename)`, etc.
- `wsprd/wsprd.h`: `const unsigned char *channel_symbols` in subtract_signal

### 4D. Clean up stale EVAL/CHECK/TODO comment markers

Remove or resolve markers at `rtlsdr_wsprd.c:174,192,729` and `wsprd/wsprd_utils.c:38`.

---

## Verification (after each phase)

```bash
make clean && make                                                         # builds OK
./rtlsdr_wsprd -f 2m -c A1XYZ -l AB12CD -t                                # self-test passes (exit 0)
./rtlsdr_wsprd -f 2m -c A1XYZ -l AB12CD -r ./signals/refSignalSnr0dB.iq   # reference IQ test passes
./rtlsdr_wsprd --version; echo $?                                          # exits 0 (after Phase 1B)
```

## Commit Strategy

1. `fix: correct self-test logic, --version exit code, and bitwise operators` (Phase 1)
2. `fix: add bounds checks, NULL guards, and volatile for thread safety` (Phase 2A-2D)
3. `fix: URL-encode user parameters in wsprnet.org reporting` (Phase 2E)
4. `fix: replace unsafe strcpy/sprintf with bounded alternatives` (Phase 2F)
5. `build: modernize Makefile with warnings, LDFLAGS, dependency tracking` (Phase 3A)
6. `ci: update outdated actions, fix typo, add .dockerignore` (Phase 3B-3C)
7. `refactor: extract constants, replace bubble sort, add const qualifiers` (Phase 4)
