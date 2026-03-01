# Plan: Add Unit and Integration Tests

## Current state

The project has only two coarse integration tests run via CLI flags:
- `./rtlsdr_wsprd -t` (synthetic self-test)
- `./rtlsdr_wsprd -r ./signals/refSignalSnr0dB.iq` (reference IQ decode)

There are no unit tests for individual functions. Many pure functions in `wsprd/` are independently testable without RTL-SDR hardware.

## Approach

No external test framework -- keep it zero-dependency. A single `tests/test_wsprd.c` file with a minimal assert-based harness (a `TEST()` macro that counts pass/fail and prints results). This matches the project's embedded/minimal style.

A shell script `tests/run_tests.sh` wraps both the unit test binary and the existing integration tests, providing a single entry point with clear pass/fail reporting.

## Test targets

### Unit tests (`tests/test_wsprd.c`)

The following pure functions can be tested with known input/output pairs:

**Message packing/unpacking round-trip** (`wsprsim_utils.c` / `wsprd_utils.c`):
- `pack_call("K1JT")` produces a known integer, `unpackcall()` recovers `"K1JT"`
- `pack_grid4_power()` / `unpackgrid()` round-trip for grid `"FN20"` + power 20
- `get_wspr_channel_symbols("K1JT FN20QI 20", ...)` produces 162 symbols, then `unpk_()` on the encoded data recovers the original callsign/grid/power

**Character code helpers** (`wsprsim_utils.c`):
- `get_callsign_character_code('A')` == 10, `('0')` == 0, `(' ')` == 36
- `get_locator_character_code('A')` == 0, `('0')` == 0, `(' ')` == 36

**Callsign unpack edge cases** (`wsprd_utils.c`):
- `unpackcall()` with value >= 262177560 returns 0 (failure)
- `unpackgrid()` with ngrid >= 32400<<7 returns 0 and writes `"XXXX"`

**Interleave/deinterleave** (`wsprsim_utils.c` / `wsprd_utils.c`):
- `interleave()` then `deinterleave()` is identity on a known 162-byte vector

**Fano encoder/decoder round-trip** (`fano.c`):
- `encode()` a known 11-byte message, then `fano()` decode recovers it exactly

**Hash function determinism** (`nhash.c`):
- `nhash("K1JT", 4, 146)` returns a stable value (compute once, hardcode as expected)

**Comparator functions** (`wsprd_utils.c`):
- `floatcomp` and `doublecomp` return correct signs for less-than, equal, greater-than

### Integration tests (`tests/run_tests.sh`)

Wraps the existing CLI tests with explicit exit-code checking:
1. Synthetic self-test: `./rtlsdr_wsprd -f 2m -c A1XYZ -l AB12CD -t`
2. Reference IQ decode: `./rtlsdr_wsprd -f 2m -c A1XYZ -l AB12CD -r ./signals/refSignalSnr0dB.iq`
3. `--version` exits 0
4. No-args prints usage and exits 0
5. Missing frequency prints error and exits 1

## Build integration

### Makefile additions

```makefile
TEST_OBJS = tests/test_wsprd.o wsprd/wsprsim_utils.o wsprd/wsprd_utils.o \
            wsprd/tab.o wsprd/fano.o wsprd/nhash.o
TARGETS = rtlsdr_wsprd tests/test_wsprd

tests/test_wsprd: $(TEST_OBJS)
	$(CC) $(LDFLAGS) -o $@ $^ -lm

.PHONY: test
test: all tests/test_wsprd
	./tests/test_wsprd
	./tests/run_tests.sh
```

The test binary links only the `wsprd/` object files plus `-lm` -- it does NOT link librtlsdr, libcurl, libfftw3, or libusb. This means tests run on any machine with just a C compiler.

### CI update (`.github/workflows/ci.yml`)

Replace the manual test commands with `make test`.

## Files to create/modify

| File | Action |
|------|--------|
| `tests/test_wsprd.c` | Create -- unit test source |
| `tests/run_tests.sh` | Create -- integration test wrapper |
| `Makefile` | Modify -- add `test` target and test binary |
| `.github/workflows/ci.yml` | Modify -- use `make test` |
| `.gitignore` | Modify -- add `tests/test_wsprd` binary |
