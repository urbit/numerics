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

// two-arg: prints `<nm> 0x<a> 0x<b> 0x<out>` (the verifier reads cols 2,3,4)
static void emit2(const char* nm, double a, double b,
                  float64_t (*fun)(float64_t, float64_t)) {
  union { double f; uint64_t b; } ua, ub; ua.f = a; ub.f = b;
  union doub xa, xb, out;
  xa.c = ua.b; xb.c = ub.b;
  out.d = fun(xa.d, xb.d);
  printf("%-6s 0x%016llx 0x%016llx 0x%016llx  (%g,%g)\n",
         nm, (unsigned long long)xa.c, (unsigned long long)xb.c,
         (unsigned long long)out.c, a, b);
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

  // trig: span several pi/2 quadrants + negatives + near-zero
  static const double tr[] = { 0.0, 0.5, 1.0, 1.5707963267948966, 2.0, 3.0,
    3.141592653589793, 4.0, 5.0, 6.283185307179586, 10.0, 100.0, 1000.0,
    -1.0, -2.0, -0.5, 0.1, -0.1, 1.0e6 };
  for (unsigned i = 0; i < sizeof tr/sizeof tr[0]; i++) emit("sin", tr[i], _rd_sin);
  for (unsigned i = 0; i < sizeof tr/sizeof tr[0]; i++) emit("cos", tr[i], _rd_cos);
  // tan: incl near-poles (close to pi/2), small, large, negatives
  static const double tn[] = { 0.0, 0.5, 1.0, 1.5, 1.5707963267948966, 2.0, 3.0,
    3.141592653589793, 4.0, 0.7, -1.0, -2.0, 0.1, -0.1, 10.0, 100.0, 0.674 };
  for (unsigned i = 0; i < sizeof tn/sizeof tn[0]; i++) emit("tan", tn[i], _rd_tan);
  // atan: spans all breakpoints (7/16,11/16,19/16,39/16) + large + negatives
  static const double at[] = { 0.0, 0.1, 0.4, 0.5, 0.6, 0.7, 1.0, 1.1, 1.5, 2.0,
    2.4, 3.0, 10.0, 100.0, 1.0e10, -0.5, -1.0, -2.0, -10.0 };
  for (unsigned i = 0; i < sizeof at/sizeof at[0]; i++) emit("atan", at[i], _rd_atan);
  // asin/acos: domain [-1,1], spanning <0.5, [0.5,0.975), near-1, +-1, tiny, OOB
  static const double iv[] = { 0.0, 0.1, 0.3, 0.5, 0.6, 0.7, 0.9, 0.975, 0.99,
    1.0, -1.0, -0.5, -0.7, -0.9, 0.25, -0.25, 2.0 };
  for (unsigned i = 0; i < sizeof iv/sizeof iv[0]; i++) emit("asin", iv[i], _rd_asin);
  for (unsigned i = 0; i < sizeof iv/sizeof iv[0]; i++) emit("acos", iv[i], _rd_acos);
  // log-2 / log-10: positive domain (reuse lg)
  for (unsigned i = 0; i < sizeof lg/sizeof lg[0]; i++) emit("log-2", lg[i], _rd_log2);
  for (unsigned i = 0; i < sizeof lg/sizeof lg[0]; i++) emit("log-10", lg[i], _rd_log10);
  // sqt: non-negative; cbt: all reals
  static const double sq[] = { 0.0, 1.0, 2.0, 4.0, 0.25, 1.0e5, 1.0e-5, 3.0,
    100.0, 1.0e300, 1.0e-300 };
  for (unsigned i = 0; i < sizeof sq/sizeof sq[0]; i++) emit("sqt", sq[i], _rd_sqt);
  static const double cb[] = { 0.0, 1.0, 2.0, 8.0, 27.0, -8.0, -2.0, 0.125,
    1.0e9, -1.0e9, 0.5, -0.5 };
  for (unsigned i = 0; i < sizeof cb/sizeof cb[0]; i++) emit("cbt", cb[i], _rd_cbt);

  // --- two-arg functions ------------------------------------------------
  // atan2(y,x): all quadrants + axes + zero cases
  static const double a2[][2] = {
    {1.0,1.0}, {1.0,-1.0}, {-1.0,1.0}, {-1.0,-1.0}, {0.0,1.0}, {0.0,-1.0},
    {1.0,0.0}, {-1.0,0.0}, {3.0,4.0}, {-3.0,4.0}, {2.0,0.5}, {0.1,100.0},
    {100.0,0.1}, {0.0,0.0} };
  for (unsigned i = 0; i < sizeof a2/sizeof a2[0]; i++)
    emit2("atan2", a2[i][0], a2[i][1], _rd_atan2);
  // pow(x,n): integer + fractional exponents, negatives, fast-path + exp/log
  static const double pw[][2] = {
    {2.0,10.0}, {2.0,0.5}, {2.0,-1.0}, {3.0,3.0}, {10.0,2.0}, {0.5,2.0},
    {2.0,3.0}, {5.0,0.0}, {9.0,0.5}, {2.0,0.1}, {1.5,4.0}, {7.0,2.0} };
  for (unsigned i = 0; i < sizeof pw/sizeof pw[0]; i++)
    emit2("pow", pw[i][0], pw[i][1], _rd_pow);
  // pow-n(x,n): positive-integer exponents (repeated mul)
  static const double pn[][2] = {
    {2.0,3.0}, {2.0,10.0}, {3.0,4.0}, {10.0,2.0}, {1.5,2.0}, {2.0,1.0},
    {5.0,3.0}, {7.0,2.0}, {0.5,4.0}, {2.0,8.0} };
  for (unsigned i = 0; i < sizeof pn/sizeof pn[0]; i++)
    emit2("pow-n", pn[i][0], pn[i][1], _rd_pow_n);

  return 0;
}
