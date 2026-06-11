// rd_check.c -- bit-exact harness for the @rd math jets.
//
// Includes the MASTER jet algorithms (jets/i/math.c) with -DMATH_JET_HARNESS,
// so the C tested here is byte-for-byte the C the runtime jet runs.  Prints
// `<fn> 0x<in> 0x<out>` per case; compare each <out> against the Hoon
// `(<fn>:rd:math <in>)` (see build.sh for the generated dojo expression).
//
// Build/run: ./build.sh   (re-archives the zig softfloat .a for Apple ld).

#define MATH_JET_HARNESS
#include "../jets/i/math.c"

#include <stdio.h>

static void emit(const char* nm, double x, float64_t (*fun)(float64_t)) {
  union { double f; uint64_t b; } u; u.f = x;
  union doub in, out;
  in.c = u.b;
  out.d = fun(in.d);
  printf("%-4s 0x%016llx 0x%016llx  (%g)\n",
         nm, (unsigned long long)in.c, (unsigned long long)out.c, x);
}

int main(void) {
  softfloat_roundingMode = softfloat_round_near_even;

  // exp: full range incl overflow/subnormal tails
  static const double ex[] = { 0.0, 0.5, 1.0, 2.0, 3.0, -1.0, -2.0, 5.0, 10.0,
    -10.0, 0.1, -0.1, 3.141592653589793, 100.0, -100.0, 700.0, -700.0, 709.0,
    -740.0 };
  for (unsigned i = 0; i < sizeof ex/sizeof ex[0]; i++) emit("exp", ex[i], _rd_exp);

  // log: positive domain, incl <1, large, small, subnormal-ish
  static const double lg[] = { 1.0, 2.0, 0.5, 10.0, 0.1, 2.718281828459045,
    100.0, 0.001, 1.0e10, 1.0e-10, 1.4142135623730951, 1.5, 3.0, 7.0,
    1.0e300, 5.0e-324 };
  for (unsigned i = 0; i < sizeof lg/sizeof lg[0]; i++) emit("log", lg[i], _rd_log);

  return 0;
}
