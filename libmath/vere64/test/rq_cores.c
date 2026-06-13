// rq_cores.c -- verify the native f128 @rq cores against math-rq.hoon vectors.
#define MATH_JET_HARNESS
#include "../noun/jets/i/math.c"
#include <stdio.h>
static int fails = 0;
static float128_t Q(c3_d hi, c3_d lo){ return _rq_bits(hi, lo); }
static void chk(const char* nm, float128_t got, c3_d ehi, c3_d elo){
  union quad u; u.q = got;
  int ok = (u.w[1]==ehi && u.w[0]==elo);
  if(!ok) fails++;
  printf("%-10s got 0x%016llx.%016llx  want 0x%016llx.%016llx  %s\n",
         nm,(unsigned long long)u.w[1],(unsigned long long)u.w[0],
         (unsigned long long)ehi,(unsigned long long)elo, ok?"OK":"*** FAIL");
}
int main(void){
  softfloat_roundingMode = softfloat_round_near_even; _math_rnd = softfloat_round_near_even;
  float128_t ONE=Q(0x3fff000000000000ULL,0), TWO=Q(0x4000000000000000ULL,0),
             HALF=Q(0x3ffe000000000000ULL,0), EIGHT=Q(0x4002000000000000ULL,0);
  chk("log-2",  _rq_log(TWO),   0x3ffe62e42fefa39eULL,0xf35793c7673007e6ULL);
  chk("log2-8", _rq_log2(EIGHT),0x4000800000000000ULL,0x0000000000000000ULL);
  chk("sin-1",  _rq_sin(ONE),   0x3ffeaed548f090ceULL,0xe0418dd3d2138a1eULL);
  chk("cos-1",  _rq_cos(ONE),   0x3ffe14a280fb5068ULL,0xb923848cdb2ed0e4ULL);
  chk("atan-1", _rq_atan(ONE),  0x3ffe921fb54442d1ULL,0x8469898cc51701b8ULL);
  chk("atan-2", _rq_atan(TWO),  0x3fff1b6e192ebbe4ULL,0x46c6d19aa220a39bULL);
  chk("asin-h", _rq_asin(HALF), 0x3ffe0c152382d736ULL,0x58465bb32e0f567bULL);
  chk("acos-h", _rq_acos(HALF), 0x3fff0c152382d736ULL,0x58465bb32e0f567bULL);
  chk("sqt-2",  _rq_sqt(TWO),   0x3fff6a09e667f3bcULL,0xc908b2fb1366ea95ULL);
  chk("cbt-8",  _rq_cbt(EIGHT), 0x3fffffffffffffffULL,0xffffffffffffffffULL);
  chk("pow-2-h",_rq_pow(TWO,HALF),0x3fff6a09e667f3bcULL,0xc908b2fb1366ea96ULL);
  chk("tan-1",  _rq_tan(ONE),   0x3fff8eb245cbee3aULL,0x5b8acc7d41323140ULL); // tan(1)=div(sin,cos) under %n
  printf("\n%s (%d failures)\n", fails?"FAILED":"ALL PASS", fails);
  return fails?1:0;
}
