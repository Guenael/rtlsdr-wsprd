/*
 * Unit tests for rtlsdr-wsprd
 *
 * Zero-dependency test harness. Links only wsprd/ objects + libm.
 * Run: ./tests/test_wsprd
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#include "../wsprd/wsprd.h"
#include "../wsprd/wsprd_utils.h"
#include "../wsprd/wsprsim_utils.h"
#include "../wsprd/fano.h"
#include "../wsprd/nhash.h"
#include "../wsprd/metric_tables.h"

static int tests_run = 0;
static int tests_passed = 0;
static int tests_failed = 0;

#define ASSERT_MSG(cond, fmt, ...)                                        \
    do {                                                                  \
        if (!(cond)) {                                                    \
            printf("  FAIL: " fmt " (%s:%d)\n", ##__VA_ARGS__,           \
                   __FILE__, __LINE__);                                   \
            return 1;                                                     \
        }                                                                 \
    } while (0)

#define ASSERT_EQ_INT(a, b) \
    ASSERT_MSG((a) == (b), "expected %d, got %d", (int)(b), (int)(a))

#define ASSERT_EQ_STR(a, b) \
    ASSERT_MSG(strcmp((a), (b)) == 0, "expected \"%s\", got \"%s\"", (b), (a))

#define ASSERT_TRUE(cond) \
    ASSERT_MSG((cond), "condition false: %s", #cond)

#define RUN_TEST(fn)                                                      \
    do {                                                                  \
        tests_run++;                                                      \
        printf("  %-50s ", #fn);                                          \
        if (fn() == 0) {                                                  \
            tests_passed++;                                               \
            printf("OK\n");                                               \
        } else {                                                          \
            tests_failed++;                                               \
        }                                                                 \
    } while (0)


/* ===== Character code helpers ===== */

static int test_callsign_character_codes(void) {
    ASSERT_EQ_INT(get_callsign_character_code('0'), 0);
    ASSERT_EQ_INT(get_callsign_character_code('9'), 9);
    ASSERT_EQ_INT(get_callsign_character_code('A'), 10);
    ASSERT_EQ_INT(get_callsign_character_code('Z'), 35);
    ASSERT_EQ_INT(get_callsign_character_code(' '), 36);
    return 0;
}

static int test_locator_character_codes(void) {
    ASSERT_EQ_INT(get_locator_character_code('0'), 0);
    ASSERT_EQ_INT(get_locator_character_code('9'), 9);
    ASSERT_EQ_INT(get_locator_character_code('A'), 0);
    ASSERT_EQ_INT(get_locator_character_code('R'), 17);
    ASSERT_EQ_INT(get_locator_character_code(' '), 36);
    return 0;
}


/* ===== pack_call / unpackcall round-trip ===== */

static int test_pack_unpack_call_k1jt(void) {
    long unsigned int n = pack_call("K1JT");
    ASSERT_TRUE(n > 0);

    char call[13] = {0};
    int ok = unpackcall((int32_t)n, call);
    ASSERT_EQ_INT(ok, 1);
    ASSERT_EQ_STR(call, "K1JT");
    return 0;
}

static int test_pack_unpack_call_va2gka(void) {
    long unsigned int n = pack_call("VA2GKA");
    ASSERT_TRUE(n > 0);

    char call[13] = {0};
    int ok = unpackcall((int32_t)n, call);
    ASSERT_EQ_INT(ok, 1);
    ASSERT_EQ_STR(call, "VA2GKA");
    return 0;
}

static int test_pack_unpack_call_w1aw(void) {
    long unsigned int n = pack_call("W1AW");
    ASSERT_TRUE(n > 0);

    char call[13] = {0};
    int ok = unpackcall((int32_t)n, call);
    ASSERT_EQ_INT(ok, 1);
    ASSERT_EQ_STR(call, "W1AW");
    return 0;
}


/* ===== unpackcall edge cases ===== */

static int test_unpackcall_out_of_range(void) {
    char call[13] = {0};
    int ok = unpackcall(262177560, call);
    ASSERT_EQ_INT(ok, 0);
    return 0;
}


/* ===== unpackgrid edge cases ===== */

static int test_unpackgrid_out_of_range(void) {
    char grid[5] = {0};
    int32_t ngrid = (32400 << 7);
    int ok = unpackgrid(ngrid, grid);
    ASSERT_EQ_INT(ok, 0);
    ASSERT_EQ_STR(grid, "XXXX");
    return 0;
}


/* ===== interleave / deinterleave identity ===== */

static int test_interleave_deinterleave_identity(void) {
    unsigned char original[162];
    unsigned char work[162];

    for (int i = 0; i < 162; i++) {
        original[i] = (unsigned char)(i & 0xFF);
        work[i] = original[i];
    }

    interleave(work);

    /* After interleave, the array should differ from original */
    int differs = 0;
    for (int i = 0; i < 162; i++) {
        if (work[i] != original[i]) differs = 1;
    }
    ASSERT_TRUE(differs);

    deinterleave(work);

    for (int i = 0; i < 162; i++) {
        ASSERT_MSG(work[i] == original[i],
                   "mismatch at index %d: expected %d, got %d",
                   i, original[i], work[i]);
    }
    return 0;
}


/* ===== Fano encode/decode round-trip ===== */

static int test_fano_encode_decode_roundtrip(void) {
    /*
     * Use the same encode path as get_wspr_channel_symbols:
     * pack a known message into 11 bytes, encode, then decode with Fano.
     * The first 7 bytes carry the 50-bit payload; bytes 7-10 are tail (zeros).
     */
    long unsigned int n = pack_call("K1JT");
    char grid4[5];
    grid4[0] = get_locator_character_code('F');
    grid4[1] = get_locator_character_code('N');
    grid4[2] = get_locator_character_code('2');
    grid4[3] = get_locator_character_code('0');
    long unsigned int m = pack_grid4_power(grid4, 20);

    unsigned char data[11] = {0};
    data[0] = 0xFF & (n >> 20);
    data[1] = 0xFF & (n >> 12);
    data[2] = 0xFF & (n >> 4);
    data[3] = ((n & 0x0F) << 4) + ((m >> 18) & 0x0F);
    data[4] = 0xFF & (m >> 10);
    data[5] = 0xFF & (m >> 2);
    data[6] = (m & 0x03) << 6;

    unsigned char encoded[11 * 8 * 2];
    memset(encoded, 0, sizeof(encoded));
    encode(encoded, data, 11);

    /* Build metric table (same as wspr_decode) */
    int32_t mettab[2][256];
    float bias = 0.45;
    for (int i = 0; i < 256; i++) {
        mettab[0][i] = roundf(10.0 * (metric_tables[2][i] - bias));
        mettab[1][i] = roundf(10.0 * (metric_tables[2][255 - i] - bias));
    }

    /* Convert hard encoder output (0/1) to soft symbols (0 or 255) */
    unsigned char soft[162];
    for (int i = 0; i < 162; i++) {
        soft[i] = encoded[i] ? 255 : 0;
    }

    unsigned char decoded[11] = {0};
    unsigned int metric, cycles, maxnp;
    int result = fano(&metric, &cycles, &maxnp, decoded, soft, 81, mettab, 60, 10000);
    ASSERT_EQ_INT(result, 0);

    for (int i = 0; i < 7; i++) {
        ASSERT_MSG(decoded[i] == data[i],
                   "byte %d: expected 0x%02X, got 0x%02X",
                   i, data[i], decoded[i]);
    }
    return 0;
}


/* ===== nhash determinism ===== */

static int test_nhash_deterministic(void) {
    uint32_t h1 = nhash("K1JT", 4, 146);
    uint32_t h2 = nhash("K1JT", 4, 146);
    ASSERT_TRUE(h1 == h2);

    /* Different inputs produce different hashes */
    uint32_t h3 = nhash("VA2GKA", 6, 146);
    ASSERT_TRUE(h1 != h3);
    return 0;
}

static int test_nhash_within_hashtab_range(void) {
    uint32_t h = nhash("K1JT", 4, 146);
    ASSERT_TRUE(h < HASHTAB_SIZE);
    return 0;
}


/* ===== Comparator functions ===== */

static int test_floatcomp(void) {
    float a = 1.0f, b = 2.0f, c = 1.0f;
    ASSERT_TRUE(floatcomp(&a, &b) < 0);
    ASSERT_TRUE(floatcomp(&b, &a) > 0);
    ASSERT_TRUE(floatcomp(&a, &c) == 0);
    return 0;
}

static int test_doublecomp(void) {
    double a = 1.0, b = 2.0, c = 1.0;
    ASSERT_TRUE(doublecomp(&a, &b) < 0);
    ASSERT_TRUE(doublecomp(&b, &a) > 0);
    ASSERT_TRUE(doublecomp(&a, &c) == 0);
    return 0;
}


/* ===== Full message encode/decode round-trip ===== */

static int test_wspr_message_roundtrip(void) {
    char message[] = "K1JT FN20QI 20";
    char hashtab[HASHTAB_SIZE * HASHTAB_ENTRY_LEN] = {0};
    char loctab[HASHTAB_SIZE * LOCTAB_ENTRY_LEN] = {0};
    unsigned char symbols[162];

    int ok = get_wspr_channel_symbols(message, hashtab, loctab, symbols);
    ASSERT_EQ_INT(ok, 1);

    /* Verify symbols are in valid range [0,3] */
    for (int i = 0; i < 162; i++) {
        ASSERT_MSG(symbols[i] <= 3,
                   "symbol[%d] = %d, expected 0-3", i, symbols[i]);
    }

    return 0;
}

static int test_wspr_message_different_inputs(void) {
    char hashtab[HASHTAB_SIZE * HASHTAB_ENTRY_LEN] = {0};
    char loctab[HASHTAB_SIZE * LOCTAB_ENTRY_LEN] = {0};
    unsigned char sym1[162], sym2[162];

    char msg1[] = "K1JT FN20QI 20";
    char msg2[] = "W1AW FN31PR 10";

    get_wspr_channel_symbols(msg1, hashtab, loctab, sym1);
    get_wspr_channel_symbols(msg2, hashtab, loctab, sym2);

    int differs = 0;
    for (int i = 0; i < 162; i++) {
        if (sym1[i] != sym2[i]) differs = 1;
    }
    ASSERT_TRUE(differs);
    return 0;
}


/* ===== pack_call edge cases ===== */

static int test_pack_call_too_long(void) {
    long unsigned int n = pack_call("TOOLONG1");
    ASSERT_EQ_INT((int)n, 0);
    return 0;
}


/* ===== unpack50 basic test ===== */

static int test_unpack50_basic(void) {
    /* Encode known n1/n2 via pack_call + pack_grid4_power, then verify unpack50 recovers them */
    long unsigned int n = pack_call("K1JT");

    char grid4[5];
    grid4[0] = get_locator_character_code('F');
    grid4[1] = get_locator_character_code('N');
    grid4[2] = get_locator_character_code('2');
    grid4[3] = get_locator_character_code('0');
    long unsigned int m = pack_grid4_power(grid4, 20);

    /* Pack into 11 bytes (same as get_wspr_channel_symbols) */
    unsigned char data[11] = {0};
    data[0] = 0xFF & (n >> 20);
    data[1] = 0xFF & (n >> 12);
    data[2] = 0xFF & (n >> 4);
    data[3] = ((n & 0x0F) << 4) + ((m >> 18) & 0x0F);
    data[4] = 0xFF & (m >> 10);
    data[5] = 0xFF & (m >> 2);
    data[6] = (m & 0x03) << 6;

    int32_t n1, n2;
    unpack50((signed char *)data, &n1, &n2);

    ASSERT_EQ_INT((int)n1, (int)n);
    ASSERT_EQ_INT((int)n2, (int)m);
    return 0;
}


/* ===== Full unpk_ round-trip ===== */

static int test_unpk_roundtrip(void) {
    char message[] = "K1JT FN20QI 20";
    char hashtab[HASHTAB_SIZE * HASHTAB_ENTRY_LEN] = {0};
    char loctab[HASHTAB_SIZE * LOCTAB_ENTRY_LEN] = {0};
    unsigned char symbols[162];

    get_wspr_channel_symbols(message, hashtab, loctab, symbols);

    /* Re-encode the packed data to get the 11-byte message back */
    long unsigned int n = pack_call("K1JT");
    char grid4[5];
    grid4[0] = get_locator_character_code('F');
    grid4[1] = get_locator_character_code('N');
    grid4[2] = get_locator_character_code('2');
    grid4[3] = get_locator_character_code('0');
    long unsigned int m = pack_grid4_power(grid4, 20);

    unsigned char data[11] = {0};
    data[0] = 0xFF & (n >> 20);
    data[1] = 0xFF & (n >> 12);
    data[2] = 0xFF & (n >> 4);
    data[3] = ((n & 0x0F) << 4) + ((m >> 18) & 0x0F);
    data[4] = 0xFF & (m >> 10);
    data[5] = 0xFF & (m >> 2);
    data[6] = (m & 0x03) << 6;

    char call_loc_pow[23] = {0};
    char call[13] = {0};
    char loc[7] = {0};
    char pwr[3] = {0};
    char callsign[13] = {0};

    int noprint = unpk_((signed char *)data, hashtab, loctab,
                        call_loc_pow, call, loc, pwr, callsign);
    ASSERT_EQ_INT(noprint, 0);
    ASSERT_EQ_STR(call, "K1JT");
    ASSERT_EQ_STR(loc, "FN20");
    ASSERT_EQ_STR(pwr, "20");
    return 0;
}


/* ===== Main ===== */

int main(void) {
    printf("=== rtlsdr-wsprd unit tests ===\n\n");

    printf("[Character codes]\n");
    RUN_TEST(test_callsign_character_codes);
    RUN_TEST(test_locator_character_codes);

    printf("\n[Pack/unpack callsign]\n");
    RUN_TEST(test_pack_unpack_call_k1jt);
    RUN_TEST(test_pack_unpack_call_va2gka);
    RUN_TEST(test_pack_unpack_call_w1aw);
    RUN_TEST(test_pack_call_too_long);
    RUN_TEST(test_unpackcall_out_of_range);

    printf("\n[Grid unpack]\n");
    RUN_TEST(test_unpackgrid_out_of_range);

    printf("\n[Unpack50]\n");
    RUN_TEST(test_unpack50_basic);

    printf("\n[Interleave/deinterleave]\n");
    RUN_TEST(test_interleave_deinterleave_identity);

    printf("\n[Fano encode/decode]\n");
    RUN_TEST(test_fano_encode_decode_roundtrip);

    printf("\n[Hash function]\n");
    RUN_TEST(test_nhash_deterministic);
    RUN_TEST(test_nhash_within_hashtab_range);

    printf("\n[Comparators]\n");
    RUN_TEST(test_floatcomp);
    RUN_TEST(test_doublecomp);

    printf("\n[WSPR message encoding]\n");
    RUN_TEST(test_wspr_message_roundtrip);
    RUN_TEST(test_wspr_message_different_inputs);

    printf("\n[Full message round-trip]\n");
    RUN_TEST(test_unpk_roundtrip);

    printf("\n--- Results: %d/%d passed, %d failed ---\n",
           tests_passed, tests_run, tests_failed);

    return tests_failed ? 1 : 0;
}
