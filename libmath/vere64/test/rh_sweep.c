// rh_sweep.c -- emit every finite f16 input->output for the native @rh kernels,
// so a Python model that mirrors the same Hoon algorithm can diff bit-for-bit
// across the WHOLE range (not just the 20 spot vectors).
#define MATH_JET_HARNESS
#include "../noun/jets/i/math.c"
#include <stdio.h>

static void sweep(const char* nm, float16_t (*fn)(float16_t)) {
  _math_rnd = softfloat_round_near_even;
  for ( unsigned b = 0; b <= 0xffff; b++ ) {
    union half u; u.c = (uint16_t)b;
    // skip NaN/inf inputs to keep the comparison about finite numerics
    unsigned exp = (b >> 10) & 0x1f, man = b & 0x3ff;
    if ( exp == 0x1f ) continue;            // inf/nan
    union half o; o.h = fn(u.h);
    printf("%s %04x %04x\n", nm, b, o.c);
  }
}
int main(void) {
  sweep("exp",  _rh_exp);
  sweep("log",  _rh_log);
  sweep("sin",  _rh_sin);
  sweep("cos",  _rh_cos);
  sweep("atan", _rh_atan);
  sweep("asin", _rh_asin);
  sweep("acos", _rh_acos);
  sweep("sqt",  _rh_sqt);
  return 0;
}
