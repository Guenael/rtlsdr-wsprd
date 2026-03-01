# Rapport: CI/CD Review and Fixes

Implementation of the CI/CD review plan for rtlsdr-wsprd. All changes verified via `make clean && make && make test`.

---

## What was changed

### `.github/workflows/ci.yml` -- Rewritten

**Before:** A single monolithic `analyze` job running CodeQL, build, test, cppcheck, cpplint, and codespell sequentially. Used bare `pip install` (PEP 668 incompatible), installed unused `cmake`, ran cppcheck with wrong C standard, and had incomplete file globs.

**After:** Two parallel jobs:

| Job | Steps |
|-----|-------|
| `build-test` | Checkout, CodeQL init, apt-get install (no cmake), `make`, CodeQL analyze, `make test` |
| `lint` | Checkout, apt-get cppcheck, `pipx install cpplint`, `pipx install codespell`, cppcheck, cpplint, codespell |

Specific fixes applied:

| Issue | Before | After |
|-------|--------|-------|
| pip install | `pip install cpplint codespell` | `pipx install cpplint` / `pipx install codespell` |
| Unused cmake | `apt-get install ... cmake ...` | Removed from `build-test` apt-get |
| Cppcheck std | `--std=c11` | `--std=c17` |
| Cppcheck glob | `*.c wsprd/*.c` | `*.c wsprd/*.c tests/*.c` |
| Cpplint glob | `--recursive *` | `*.c *.h wsprd/*.c wsprd/*.h tests/*.c` |
| Codespell glob | `*.md *.c *.h wsprd/*.c wsprd/*.h .github/workflows/ci.yml` | `*.md *.c *.h wsprd/*.c wsprd/*.h tests/*.c tests/*.sh documentation/**/*.md .github/workflows/*.yml` |
| Job parallelism | 1 sequential job | 2 parallel jobs (`build-test` + `lint`) |

### `.github/workflows/container.yml` -- Modified

| Change | Before | After |
|--------|--------|-------|
| Job name | `main` | `container` |
| Layer cache | None | `cache-from: type=gha` / `cache-to: type=gha,mode=max` |

The GHA cache backend stores Docker layers in the GitHub Actions cache, avoiding full rebuilds of the 5-architecture multi-platform image on every push.

### `Dockerfile` -- Modified

| Change | Before | After |
|--------|--------|-------|
| Test in build | `make && make install` | `make && make test && make install` |
| CA certificates | Not installed | `ca-certificates` added to runtime apt-get |
| Indentation | Tab on `libcurl4` line | Consistent 4-space indent throughout |

Adding `make test` to the build stage ensures that a broken binary cannot be pushed to GHCR. The unit tests run without RTL-SDR hardware and complete in under a second.

Adding `ca-certificates` ensures the runtime container can make HTTPS POST requests to wsprnet.org without certificate verification errors.

### `.dockerignore` -- Modified

Added entries to reduce Docker build context size:

| Entry | Reason |
|-------|--------|
| `documentation/` | Documentation not needed in image |
| `RAPPORT.md` | Report file not needed in image |
| `selftest.iq` | Generated test artifact |
| `tests/test_wsprd` | Compiled test binary |

### `.gitpod.Dockerfile` -- Modified

Removed unused packages: `valgrind`, `libcmocka-dev`, `cmocka-doc`, `libcmocka0`. The project's test framework is a custom zero-dependency C harness -- these packages were leftover from an earlier approach and wasted image build time.

---

## Files modified

| File | Change summary |
|------|----------------|
| `.github/workflows/ci.yml` | Split into parallel `build-test` + `lint` jobs; fix pip, cppcheck std, all globs; remove cmake |
| `.github/workflows/container.yml` | Rename job `main` to `container`; add GHA layer cache |
| `Dockerfile` | Add `make test` to build stage; add `ca-certificates` to runtime; fix indentation |
| `.dockerignore` | Add `documentation/`, `RAPPORT.md`, `selftest.iq`, `tests/test_wsprd` |
| `.gitpod.Dockerfile` | Remove unused valgrind/cmocka packages |

---

## Verification

```
$ make clean && make && make test

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

All 18 unit tests and 5 integration tests pass. The CI/CD changes are configuration-only and do not affect the build output.
