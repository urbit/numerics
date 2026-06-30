// rd_check.c -- bit-exact harness for the @rd math jets.
//
// Includes the MASTER jet algorithms (noun/jets/i/math.c) with
// -DMATH_JET_HARNESS, so the C tested here is byte-for-byte the C the runtime
// jet runs.  Prints `<fn> 0x<in> 0x<out>` per case; compare each <out> against
// the Hoon `(<fn>:rd:math <in>)` (see build.sh for the generated dojo
// expression).  The shared `_rd_*` cores are word-size-agnostic, so this 64-bit
// copy and the 32-bit libmath/vere/ copy produce identical results.
//
// Build/run: ./build.sh   (re-archives the zig softfloat .a for Apple ld).

#define MATH_JET_HARNESS
#include "../noun/jets/i/math.c"

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

// @rs single-arg: input given as the 32-bit bit pattern directly.
static void emits(const char* nm, unsigned in, float32_t (*fun)(float32_t)) {
  union sing x, out;
  x.c = (uint32_t)in;
  out.s = fun(x.s);
  printf("%-6s 0x%08x 0x%08x\n", nm, (unsigned)x.c, (unsigned)out.c);
}

// narrow spike: f32 -> f16 via SoftFloat in each rounding mode; tag is the
// mode letter so the Hoon (~(narrow-sh rh [<mode> ..]) <f32>) can be compared.
static void emit_narrow(unsigned f32bits) {
  union sing x; x.c = (uint32_t)f32bits;
  const struct { char m; int sf; } modes[] = {
    {'n', softfloat_round_near_even}, {'z', softfloat_round_minMag},
    {'u', softfloat_round_max}, {'d', softfloat_round_min} };
  for (unsigned k = 0; k < 4; k++) {
    softfloat_roundingMode = modes[k].sf;
    float16_t h = f32_to_f16(x.s);
    union { float16_t h; uint16_t c; } o; o.h = h;
    printf("narrow-%c 0x%08x 0x%04x\n", modes[k].m, (unsigned)x.c, (unsigned)o.c);
  }
  softfloat_roundingMode = softfloat_round_near_even;
}

int main(void) {
  softfloat_roundingMode = softfloat_round_near_even;
  // --- @rq exp probe: print hi.lo of _rq_exp(hi,lo) ---
  { struct { c3_d hi, lo; } qx[] = {
      {0x3fff000000000000ULL,0}, {0x4000000000000000ULL,0}, {0,0},
      {0xbfff000000000000ULL,0}, {0x4002000000000000ULL,0},        // -1, ~ exp args
      {0x3ffb999999999999ULL,0x999999999999999aULL} };             // 0.1
    for (unsigned i=0;i<sizeof qx/sizeof qx[0];i++){
      union quad in,out; in.w[1]=qx[i].hi; in.w[0]=qx[i].lo;
      out.q=_rq_exp(in.q);
      printf("exp-q 0x%016llx.%016llx 0x%016llx.%016llx\n",
        (unsigned long long)in.w[1],(unsigned long long)in.w[0],
        (unsigned long long)out.w[1],(unsigned long long)out.w[0]);
    } }
  // --- @rh narrow spike: stress f32->f16 in all 4 modes ---
  static const unsigned nrw[] = {
    0x3f800000, 0x40490fdb, 0xc0490fdb, 0x477fe000, 0x477ff000, 0x477ff001,
    0x47800000, 0x7f7fffff, 0xff7fffff, 0x33800000, 0x33000000, 0x33000001,
    0x387fc000, 0x38800000, 0x38801000, 0x33c00000, 0x7f800000, 0xff800000,
    0x7fc00000, 0x00000000, 0x80000000, 0x38ffe000, 0x477fefff };
  for (unsigned i = 0; i < sizeof nrw/sizeof nrw[0]; i++) emit_narrow(nrw[i]);
  // @rs exp: core + edges (expected from math.hoon tests/lib/math-exp.hoon)
  static const unsigned rsex[] = { 0x0, 0x3f000000, 0x3f800000, 0xbf800000,
    0x40000000, 0x41200000, 0xc0a00000, 0x3dcccccd, 0x7f800000, 0xff800000,
    0x7fc00000, 0x42b00000, 0x42b20000, 0x42c80000, 0xc2ce0000, 0xc2d00000,
    0xc2dc0000 };
  for (unsigned i = 0; i < sizeof rsex/sizeof rsex[0]; i++)
    emits("exp-s", rsex[i], _rs_exp);
  static const unsigned rstr[] = { 0x0, 0x3f000000, 0x3f800000, 0xbf800000,
    0x41200000, 0x42c80000, 0x7f800000, 0xff800000, 0x80000000, 0x40000000,
    0x40490fdb, 0x42652ee0 };
  for (unsigned i = 0; i < sizeof rstr/sizeof rstr[0]; i++) emits("sin-s", rstr[i], _rs_sin);
  for (unsigned i = 0; i < sizeof rstr/sizeof rstr[0]; i++) emits("cos-s", rstr[i], _rs_cos);
  for (unsigned i = 0; i < sizeof rstr/sizeof rstr[0]; i++) emits("tan-s", rstr[i], _rs_tan);
  // atan / asin / acos / log / log-2 / log-10 / sqt / cbt (expected: math-*.hoon)
  static const unsigned at_s[] = { 0x3f000000, 0x3f800000, 0xbf800000, 0x40000000,
    0x41200000, 0x7f800000, 0x7fc00000, 0xff800000, 0x0, 0x40490fdb, 0x3dcccccd };
  for (unsigned i=0;i<sizeof at_s/sizeof at_s[0];i++) emits("atan-s", at_s[i], _rs_atan);
  static const unsigned iv_s[] = { 0x0, 0x3f000000, 0x3f800000, 0xbf800000, 0x3f400000,
    0x3f666666, 0xbf19999a, 0x7fc00000, 0x3dcccccd };
  for (unsigned i=0;i<sizeof iv_s/sizeof iv_s[0];i++) emits("asin-s", iv_s[i], _rs_asin);
  for (unsigned i=0;i<sizeof iv_s/sizeof iv_s[0];i++) emits("acos-s", iv_s[i], _rs_acos);
  static const unsigned lg_s[] = { 0x3f800000, 0x40000000, 0x3f000000, 0x41200000,
    0x42c80000, 0x3dcccccd, 0x40ec7326, 0x000116c2, 0x7f800000, 0xff800000,
    0x7fc00000, 0x0, 0xbf800000, 0x41000000, 0x447a0000, 0x40c90fdb };
  for (unsigned i=0;i<sizeof lg_s/sizeof lg_s[0];i++) emits("log-s", lg_s[i], _rs_log);
  for (unsigned i=0;i<sizeof lg_s/sizeof lg_s[0];i++) emits("log2-s", lg_s[i], _rs_log2);
  for (unsigned i=0;i<sizeof lg_s/sizeof lg_s[0];i++) emits("log10-s", lg_s[i], _rs_log10);
  static const unsigned sq_s[] = { 0x40000000, 0x3f000000, 0x40800000, 0x47c35000,
    0x3dcccccd, 0x41100000, 0x0, 0x7f800000, 0xbf800000 };
  for (unsigned i=0;i<sizeof sq_s/sizeof sq_s[0];i++) emits("sqt-s", sq_s[i], _rs_sqt);
  static const unsigned cb_s[] = { 0x41000000, 0xc1000000, 0x40000000, 0x0, 0x3f800000,
    0xbf800000, 0x40400000 };
  for (unsigned i=0;i<sizeof cb_s/sizeof cb_s[0];i++) emits("cbt-s", cb_s[i], _rs_cbt);

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
    {2.0,10.0}, {2.0,0.5}, {2.0,-1.0}, {3.0,3.0}, {10.0,2.0}, {0.5,2.0}, {3.0,2.5},
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
