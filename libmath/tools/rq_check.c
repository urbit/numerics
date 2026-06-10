// rq_check.c -- SoftFloat-f128 reference for @rq transcendentals (PR #18).
//
// SoftFloat = the jet arithmetic (this code IS the jet body, modulo the u3 ABI);
// MPFR = correctly-rounded truth.  Coeffs from mpmath chebyfit @113-bit.
//
// Build (macOS/arm64; SoftFloat is the f128M pointer API, vendored in vere):
//   INC=$(dirname $(find ~/urbit/vere/zig-pkg -name softfloat.h | head -1))
//   SF=$(find ~/urbit/vere -name libsoftfloat.a | head -1)
//   # the zig .a is not 8-byte aligned for Apple ld; re-archive it:
//   mkdir -p /tmp/sfo && (cd /tmp/sfo && ar x "$SF" && chmod u+rw *.o)
//   libtool -static -o /tmp/libsoftfloat.a /tmp/sfo/*.o
//   cc rq_check.c -I"$INC" -I/opt/homebrew/include /tmp/libsoftfloat.a \
//      -L/opt/homebrew/lib -lmpfr -lgmp -o /tmp/rq_check && /tmp/rq_check
//
// exp @rq: 0.998 ULP over [-20,20] (MPFR truth); output matches cheb_check.py
// (the mpmath f128 model) bit-for-bit.
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include "softfloat.h"
#include <mpfr.h>
typedef float128_t q;
static q qadd(q a,q b){q r; f128M_add(&a,&b,&r); return r;}
static q qsub(q a,q b){q r; f128M_sub(&a,&b,&r); return r;}
static q qmul(q a,q b){q r; f128M_mul(&a,&b,&r); return r;}
static q qi(int64_t n){q r; i64_to_f128M(n,&r); return r;}
static q q2k(int k){q r; r.v[0]=0; r.v[1]=((uint64_t)(k+16383))<<48; return r;}
static void qprint(const char*nm,q a){printf("%s0x%016llx%016llx\n",nm,(unsigned long long)a.v[1],(unsigned long long)a.v[0]);}
static void q_to_mpfr(mpfr_t o,q a){uint64_t hi=a.v[1],lo=a.v[0]; int s=hi>>63; int e=(hi>>48)&0x7fff;
  uint64_t mh=hi&0xffffffffffffULL; mpfr_set_ui(o,mh,MPFR_RNDN); mpfr_mul_2ui(o,o,64,MPFR_RNDN);
  mpfr_t t; mpfr_init2(t,128); mpfr_set_ui(t,lo,MPFR_RNDN); mpfr_add(o,o,t,MPFR_RNDN); mpfr_clear(t);
  mpfr_t one; mpfr_init2(one,128); mpfr_set_ui(one,1,MPFR_RNDN); mpfr_mul_2ui(one,one,112,MPFR_RNDN);
  mpfr_add(o,o,one,MPFR_RNDN); mpfr_clear(one); mpfr_mul_2si(o,o,e-16383-112,MPFR_RNDN); if(s)mpfr_neg(o,o,MPFR_RNDN);}
static const q LOG2E={{0xe1777d0ffda0d23aull,0x3fff71547652b82full}}, LN2HI={{0xf35793c600000000ull,0x3ffe62e42fefa39eull}}, LN2LO={{0x81e6864ce5316c5bull,0x3fae673007e5ed5eull}}, HALF={{0x0000000000000000ull,0x3ffe000000000000ull}};
static const q EXC[25]={
  {{0x0000000000000000ull,0x3fff000000000000ull}}, {{0x0000000000000000ull,0x3fff000000000000ull}}, {{0x0000000000000000ull,0x3ffe000000000000ull}}, {{0x5555555555555555ull,0x3ffc555555555555ull}}, {{0x5555555555555555ull,0x3ffa555555555555ull}}, {{0x1111111111111111ull,0x3ff8111111111111ull}}, {{0x6c16c16c16c16c17ull,0x3ff56c16c16c16c1ull}}, {{0xa01a01a01a01a3e8ull,0x3ff2a01a01a01a01ull}}, {{0xa01a01a01a01a146ull,0x3fefa01a01a01a01ull}}, {{0x38faac1c88a5a526ull,0x3fec71de3a556c73ull}}, {{0xc72ef016d3d6e867ull,0x3fe927e4fb7789f5ull}}, {{0x38fe748363c46e8bull,0x3fe5ae64567f544eull}}, {{0x7b544dab18f475c5ull,0x3fe21eed8eff8d89ull}}, {{0x97c9f3aebabb2423ull,0x3fde6124613a86d0ull}}, {{0xd20b83c7f94d17d8ull,0x3fda93974a8c07c9ull}}, {{0xf5f4284f0d74f9e7ull,0x3fd6ae7f3e733b81ull}}, {{0xf417b4d27c5f92a9ull,0x3fd2ae7f3e733b81ull}}, {{0x6a419e674779c97cull,0x3fce952c77030a99ull}}, {{0x0466ff8c8b42b3dfull,0x3fca6827863b97b5ull}}, {{0x874b7a686d819241ull,0x3fc62f49b469f892ull}}, {{0xbb3b32a11bb5f139ull,0x3fc1e542ba427463ull}}, {{0xc93890ff9ab55cbbull,0x3fbd71b8db9f7f73ull}}, {{0x6efc0717eae785a1ull,0x3fb90ce38aab7bd7ull}}, {{0xcb3f4f7edfaa2666ull,0x3fb47693274bab2aull}}, {{0x61cb0e23655d47cbull,0x3faff3629154e0a7ull}}
};
static q expq(q x){
  q xl=qmul(x,LOG2E), t=qadd(xl,HALF);
  int64_t k=f128M_to_i64(&t,softfloat_round_min,false); q qk=qi(k);
  q r=qsub(qsub(x,qmul(qk,LN2HI)),qmul(qk,LN2LO));
  q p=EXC[24]; for(int i=24-1;i>=0;i--) p=qadd(qmul(p,r),EXC[i]);
  return qmul(p,q2k((int)k));
}
static double ulp_err(q g,mpfr_t tr){
  mpfr_t gv,d,u; mpfr_inits2(300,gv,d,u,(mpfr_ptr)0); q_to_mpfr(gv,g);
  mpfr_sub(d,gv,tr,MPFR_RNDN); long e2; mpfr_get_d_2exp(&e2,gv,MPFR_RNDN);
  mpfr_set_ui(u,1,MPFR_RNDN); mpfr_mul_2si(u,u,(int)e2-1-112,MPFR_RNDN); mpfr_div(d,d,u,MPFR_RNDN);
  double r=mpfr_get_d(d,MPFR_RNDN); if(r<0)r=-r; mpfr_clears(gv,d,u,(mpfr_ptr)0); return r;}
int main(void){
  softfloat_roundingMode=softfloat_round_near_even;
  mpfr_t tr; mpfr_init2(tr,300); double worst=0;
  for(int i=-2000;i<=2000;i++){ q x=qmul(qi(i),qmul(qi(1),HALF)); // i*0.5 in [-1000,1000]? use i/100
    (void)x; }
  for(int i=-2000;i<=2000;i++){ q x=qsub(qmul(qi(i),HALF),qmul(qi(i),qsub(HALF,qmul(HALF,qmul(HALF,qi(1)))))); (void)x; }
  // simple sweep: x = i/128 for i in [-2560,2560] -> [-20,20]
  q inv128=q2k(-7);
  for(int i=-2560;i<=2560;i++){ q x=qmul(qi(i),inv128); q g=expq(x);
    mpfr_set_si(tr,i,MPFR_RNDN); mpfr_div_ui(tr,tr,128,MPFR_RNDN); mpfr_exp(tr,tr,MPFR_RNDN);
    double e=ulp_err(g,tr); if(e>worst)worst=e; }
  printf("# rq exp via SoftFloat f128 (jet basis); max %.3f ULP over [-20,20]\n",worst);
  qprint("exp(1)  = ",expq(qi(1)));
  qprint("exp(.5) = ",expq(HALF));
  qprint("exp(-2) = ",expq(qi(-2)));
  qprint("exp(10) = ",expq(qi(10)));
  mpfr_clear(tr); return 0;
}
