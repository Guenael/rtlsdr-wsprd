# Plan: CI/CD Review and Fixes

## Current state

The CI/CD pipeline has accumulated several issues ranging from deprecated tooling to missing safety nets. A full audit of `.github/workflows/ci.yml`, `.github/workflows/container.yml`, `Dockerfile`, `.dockerignore`, and `.gitpod.Dockerfile` was performed.

## Findings

### ci.yml (7 issues)

| # | Severity | Issue |
|---|----------|-------|
| 1 | medium | Bare `pip install cpplint codespell` -- blocked by PEP 668 on Ubuntu 24.04+ runners |
| 2 | low | `cmake` installed but never used (project uses Makefiles) |
| 3 | medium | `cppcheck --std=c11` mismatches the Makefile's `-std=gnu17` |
| 4 | low | Cppcheck glob misses `tests/*.c` |
| 5 | low | Codespell glob misses `tests/*.c`, `tests/*.sh`, `documentation/**/*.md` |
| 6 | low | `cpplint --recursive *` scans non-source directories |
| 7 | medium | Single monolithic job -- lint waits for build+CodeQL to finish |

### container.yml (2 issues)

| # | Severity | Issue |
|---|----------|-------|
| 8 | low | Job name `main` is confusing (same as the branch name) |
| 9 | medium | No Docker layer cache -- 5-architecture build from scratch every time |

### Dockerfile (3 issues)

| # | Severity | Issue |
|---|----------|-------|
| 10 | medium | No `make test` in build stage -- broken binary could be pushed to GHCR |
| 11 | medium | Missing `ca-certificates` in runtime image -- HTTPS to wsprnet.org fails |
| 12 | cosmetic | Mixed tab/space indentation in runtime `apt-get` block |

### .dockerignore (1 issue)

| # | Severity | Issue |
|---|----------|-------|
| 13 | low | Missing `documentation/`, `RAPPORT.md`, `selftest.iq`, `tests/test_wsprd` -- bloats build context |

### .gitpod.Dockerfile (1 issue)

| # | Severity | Issue |
|---|----------|-------|
| 14 | low | Installs unused `valgrind`, `libcmocka-dev`, `cmocka-doc`, `libcmocka0` |

## Proposed changes

### ci.yml

- Split into two parallel jobs: `build-test` (build + CodeQL + `make test`) and `lint` (cppcheck, cpplint, codespell). They run concurrently, cutting wall-clock time roughly in half.
- Replace bare `pip install` with `pipx install` (PEP 668 safe).
- Remove `cmake` from apt-get.
- Fix cppcheck `--std=c11` to `--std=c17`, add `tests/*.c`.
- Expand codespell globs to cover `tests/*.c`, `tests/*.sh`, `documentation/**/*.md`, `.github/workflows/*.yml`.
- Target cpplint at explicit source globs: `*.c *.h wsprd/*.c wsprd/*.h tests/*.c`.

### container.yml

- Rename job from `main` to `container`.
- Add Docker layer caching via `cache-from: type=gha` / `cache-to: type=gha,mode=max`.

### Dockerfile

- Add `make test` to the build stage (unit tests don't need hardware).
- Add `ca-certificates` to runtime stage for HTTPS.
- Fix mixed tab/space indentation.

### .dockerignore

- Add `documentation/`, `RAPPORT.md`, `selftest.iq`, `tests/test_wsprd`.

### .gitpod.Dockerfile

- Remove unused `valgrind`, `libcmocka-dev`, `cmocka-doc`, `libcmocka0`.

## Files to modify

| File | Action |
|------|--------|
| `.github/workflows/ci.yml` | Rewrite -- split jobs, fix pip, cppcheck, globs |
| `.github/workflows/container.yml` | Modify -- rename job, add cache |
| `Dockerfile` | Modify -- add test + ca-certificates, fix indent |
| `.dockerignore` | Modify -- add exclusions |
| `.gitpod.Dockerfile` | Modify -- remove unused packages |
