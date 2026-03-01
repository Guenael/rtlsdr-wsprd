# Rapport: Unit and Integration Tests

Implementation of the test plan for rtlsdr-wsprd. All tests verified passing via `make test`.

---

## What was added

### Unit tests (`tests/test_wsprd.c`)

A zero-dependency C test harness with 18 tests covering the pure-function layers of the WSPR decoder. The test binary links only `wsprd/` objects + `-lm` -- no RTL-SDR, curl, FFTW, or USB libraries required.

| # | Test | What it verifies |
|---|------|------------------|
| 1 | `test_callsign_character_codes` | `get_callsign_character_code()` maps '0'->0, '9'->9, 'A'->10, 'Z'->35, ' '->36 |
| 2 | `test_locator_character_codes` | `get_locator_character_code()` maps '0'->0, '9'->9, 'A'->0, 'R'->17, ' '->36 |
| 3 | `test_pack_unpack_call_k1jt` | `pack_call("K1JT")` round-trips through `unpackcall()` |
| 4 | `test_pack_unpack_call_va2gka` | `pack_call("VA2GKA")` round-trips through `unpackcall()` |
| 5 | `test_pack_unpack_call_w1aw` | `pack_call("W1AW")` round-trips through `unpackcall()` |
| 6 | `test_pack_call_too_long` | `pack_call("TOOLONG1")` returns 0 (>6 chars rejected) |
| 7 | `test_unpackcall_out_of_range` | `unpackcall(262177560)` returns 0 (boundary rejection) |
| 8 | `test_unpackgrid_out_of_range` | `unpackgrid(32400<<7)` returns 0 and writes "XXXX" |
| 9 | `test_unpack50_basic` | Pack known K1JT/FN20/20 into 11 bytes, verify `unpack50()` recovers n1 and n2 |
| 10 | `test_interleave_deinterleave_identity` | `interleave()` then `deinterleave()` is identity on 162-byte vector |
| 11 | `test_fano_encode_decode_roundtrip` | `encode()` a WSPR message, `fano()` decode recovers all 7 payload bytes |
| 12 | `test_nhash_deterministic` | `nhash("K1JT", 4, 146)` returns same value on repeated calls; different input gives different hash |
| 13 | `test_nhash_within_hashtab_range` | `nhash()` output is < `HASHTAB_SIZE` (32768) |
| 14 | `test_floatcomp` | `floatcomp()` returns correct signs for <, >, == |
| 15 | `test_doublecomp` | `doublecomp()` returns correct signs for <, >, == |
| 16 | `test_wspr_message_roundtrip` | `get_wspr_channel_symbols("K1JT FN20QI 20", ...)` produces 162 symbols all in [0,3] |
| 17 | `test_wspr_message_different_inputs` | Two different messages produce different symbol vectors |
| 18 | `test_unpk_roundtrip` | Full encode-then-`unpk_()` round-trip recovers call="K1JT", loc="FN20", pwr="20" |

### Integration tests (`tests/run_tests.sh`)

A shell script wrapping 5 CLI-level tests with explicit exit-code checking:

| # | Test | What it verifies |
|---|------|------------------|
| 1 | `--version exits 0` | `--version` prints version and exits successfully |
| 2 | `no args prints usage, exits 0` | Running with no arguments shows usage and exits 0 |
| 3 | `missing callsign exits non-zero` | `-f 2m` alone (no `-c`) exits with failure |
| 4 | `synthetic self-test (-t)` | Full encode/decode self-test passes |
| 5 | `decode refSignalSnr0dB.iq` | Reference IQ file decodes successfully |

### Build system changes

**`Makefile`:**
- Added `TEST_OBJS` variable for the test binary's object files
- Added `tests/test_wsprd` build target linking only `wsprd/` objects + `-lm`
- Added `test` phony target that builds everything then runs both test suites
- Updated `clean` to remove `tests/*.o`, `tests/*.d`, and `tests/test_wsprd`

**`.github/workflows/ci.yml`:**
- Replaced manual test commands with `make test`

**`.gitignore`:**
- Added `tests/test_wsprd` binary

---

## Files created

| File | Description |
|------|-------------|
| `tests/test_wsprd.c` | Unit test source (18 tests, zero dependencies beyond wsprd/ + libm) |
| `tests/run_tests.sh` | Integration test shell script (5 tests) |

## Files modified

| File | Change |
|------|--------|
| `Makefile` | Added `TEST_OBJS`, test binary target, `test` phony target, updated `clean` |
| `.github/workflows/ci.yml` | Replaced manual test commands with `make test` |
| `.gitignore` | Added `tests/test_wsprd` |

---

## Verification

```
$ make clean && make CC=gcc test

./tests/test_wsprd
=== rtlsdr-wsprd unit tests ===

[Character codes]
  test_callsign_character_codes                      OK
  test_locator_character_codes                       OK

[Pack/unpack callsign]
  test_pack_unpack_call_k1jt                         OK
  test_pack_unpack_call_va2gka                       OK
  test_pack_unpack_call_w1aw                         OK
  test_pack_call_too_long                            OK
  test_unpackcall_out_of_range                       OK

[Grid unpack]
  test_unpackgrid_out_of_range                       OK

[Unpack50]
  test_unpack50_basic                                OK

[Interleave/deinterleave]
  test_interleave_deinterleave_identity              OK

[Fano encode/decode]
  test_fano_encode_decode_roundtrip                  OK

[Hash function]
  test_nhash_deterministic                           OK
  test_nhash_within_hashtab_range                    OK

[Comparators]
  test_floatcomp                                     OK
  test_doublecomp                                    OK

[WSPR message encoding]
  test_wspr_message_roundtrip                        OK
  test_wspr_message_different_inputs                 OK

[Full message round-trip]
  test_unpk_roundtrip                                OK

--- Results: 18/18 passed, 0 failed ---

./tests/run_tests.sh
=== rtlsdr-wsprd integration tests ===

[CLI basics]
  --version exits 0                                  OK
  no args prints usage, exits 0                      OK
  missing callsign exits non-zero                    OK

[Decoder self-test]
  synthetic self-test (-t)                           OK

[Reference IQ decode]
  decode refSignalSnr0dB.iq                          OK

--- Results: 5 passed, 0 failed ---
```

## Design notes

- **No external dependencies.** The test harness uses a simple `RUN_TEST` / `ASSERT_*` macro pattern. No CMake, no CTest, no Unity, no Check. Just C and a shell script.
- **Fast.** The unit test binary runs in under 50ms. The full `make test` including integration tests completes in about 1 second.
- **Minimal link footprint.** The test binary links only the `wsprd/` codec objects and `-lm`. It does not pull in librtlsdr, libcurl, libfftw3, or libusb, so it builds and runs on any machine with a C compiler.
- **CI-integrated.** `make test` is the single entry point used by both developers and GitHub Actions.
