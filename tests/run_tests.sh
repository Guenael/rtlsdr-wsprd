#!/bin/bash
#
# Integration tests for rtlsdr_wsprd
# Requires the binary to be built first (make all).
#

set -e

BINARY="./rtlsdr_wsprd"
PASS=0
FAIL=0

run_expect_success() {
    local desc="$1"
    shift
    printf "  %-50s " "$desc"
    if "$@" > /dev/null 2>&1; then
        echo "OK"
        PASS=$((PASS + 1))
    else
        echo "FAIL (exit $?)"
        FAIL=$((FAIL + 1))
    fi
}

run_expect_failure() {
    local desc="$1"
    shift
    printf "  %-50s " "$desc"
    if "$@" > /dev/null 2>&1; then
        echo "FAIL (expected non-zero exit)"
        FAIL=$((FAIL + 1))
    else
        echo "OK"
        PASS=$((PASS + 1))
    fi
}

echo "=== rtlsdr-wsprd integration tests ==="
echo

if [ ! -x "$BINARY" ]; then
    echo "ERROR: $BINARY not found or not executable. Run 'make' first."
    exit 1
fi

echo "[CLI basics]"
run_expect_success "--version exits 0" \
    $BINARY --version

run_expect_success "no args prints usage, exits 0" \
    $BINARY

run_expect_failure "missing callsign exits non-zero" \
    $BINARY -f 2m

echo
echo "[Decoder self-test]"
run_expect_success "synthetic self-test (-t)" \
    $BINARY -f 2m -c A1XYZ -l AB12CD -t

echo
echo "[Reference IQ decode]"
run_expect_success "decode refSignalSnr0dB.iq" \
    $BINARY -f 2m -c A1XYZ -l AB12CD -r ./signals/refSignalSnr0dB.iq

echo
echo "--- Results: $PASS passed, $FAIL failed ---"

[ "$FAIL" -eq 0 ] || exit 1
