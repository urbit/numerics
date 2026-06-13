// rh_check.c -- verify the native f16 @rh cores against the math-rh.hoon test
// vectors.  Includes the master math.c with -DMATH_JET_HARNESS so the C tested
// is byte-for-byte the C the runtime compiles.  Door r=%n in the tests, so we
// set _math_rnd = near_even before every call.
#define MATH_JET_HARNESS
#include "../noun/jets/i/math.c"

#include <stdio.h>

static int fails = 0;
// Model the real jet wrapper: pre-dirty the global rounding mode to %z (as a
// prior door op would leave it on-ship), then set near-even exactly as the
// fixed _rh_jet/u3qi wrappers do.  This catches the "wrapper set _math_rnd but
// not softfloat_roundingMode" regression.
static uint16_t b1(uint16_t in, float16_t (*fn)(float16_t)) {
  softfloat_roundingMode = softfloat_round_minMag;     // simulate stale %z
  softfloat_roundingMode = softfloat_round_near_even;  // wrapper resets it
  _math_rnd = softfloat_round_near_even;
  union half u; u.h = fn(_rh_bits(in).h); return u.c;
}
static uint16_t b2(uint16_t a, uint16_t b, float16_t (*fn)(float16_t, float16_t)) {
  softfloat_roundingMode = softfloat_round_minMag;
  softfloat_roundingMode = softfloat_round_near_even;
  _math_rnd = softfloat_round_near_even;
  union half u; u.h = fn(_rh_bits(a).h, _rh_bits(b).h); return u.c;
}
static void chk(const char* nm, uint16_t got, uint16_t want) {
  int ok = (got == want);
  if ( !ok ) fails++;
  printf("%-12s got 0x%04x  want 0x%04x  %s\n", nm, got, want, ok ? "OK" : "*** FAIL");
}

int main(void) {
  chk("exp-half",  b1(0x3800, _rh_exp),  0x3e98);
  chk("exp-1",     b1(0x3c00, _rh_exp),  0x4170);
  chk("exp-n2",    b1(0xc000, _rh_exp),  0x3055);
  chk("exp-inf",   b1(0x7c00, _rh_exp),  0x7c00);
  chk("log-2",     b1(0x4000, _rh_log),  0x398c);
  chk("log-half",  b1(0x3800, _rh_log),  0xb98c);
  chk("sin-1",     b1(0x3c00, _rh_sin),  0x3abb);
  chk("sin-pi",    b1(0x4248, _rh_sin),  0x13ed);
  chk("cos-1",     b1(0x3c00, _rh_cos),  0x3853);
  chk("tan-1",     b1(0x3c00, _rh_tan),  0x3e3a);
  chk("atan-1",    b1(0x3c00, _rh_atan), 0x3a48);
  chk("atan-2",    b1(0x4000, _rh_atan), 0x3c6e);
  chk("asin-half", b1(0x3800, _rh_asin), 0x3830);
  chk("acos-half", b1(0x3800, _rh_acos), 0x3c30);
  chk("sqt-2",     b1(0x4000, _rh_sqt),  0x3da8);
  chk("sqt-10",    b1(0x4900, _rh_sqt),  0x4253);
  chk("cbt-8",     b1(0x4800, _rh_cbt),  0x4000);
  chk("log2-8",    b1(0x4800, _rh_log2), 0x4200);
  chk("log10-1k",  b1(0x6400, _rh_log10),0x4205);
  chk("pow-2-h",   b2(0x4000, 0x3800, _rh_pow), 0x3da8);
  printf("\n%s (%d failures)\n", fails ? "FAILED" : "ALL PASS", fails);
  return fails ? 1 : 0;
}
