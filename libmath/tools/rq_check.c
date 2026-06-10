// rq_check.c -- SoftFloat-f128 reference for @rq transcendentals (PR #18).
// SoftFloat = the jet arithmetic (this code IS the jet body, modulo the u3 ABI);
// MPFR = correctly-rounded truth.  Coeffs from mpmath chebyfit/Taylor @113-bit.
// Build recipe: re-archive libsoftfloat.a with `libtool -static` (zig .a is not
// 8-byte aligned for Apple ld), then cc -I<sf-include> -I/opt/homebrew/include
// rq_check.c /tmp/libsoftfloat.a -lmpfr -lgmp.
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include "softfloat.h"
#include <mpfr.h>
typedef float128_t q;
static q qadd(q a,q b){q r; f128M_add(&a,&b,&r); return r;}
static q qsub(q a,q b){q r; f128M_sub(&a,&b,&r); return r;}
static q qmul(q a,q b){q r; f128M_mul(&a,&b,&r); return r;}
static q qdiv(q a,q b){q r; f128M_div(&a,&b,&r); return r;}
static q qi(int64_t n){q r; i64_to_f128M(n,&r); return r;}
static q q2k(int k){q r; r.v[0]=0; r.v[1]=((uint64_t)(k+16383))<<48; return r;}
static void qprint(const char*nm,q a){printf("%s0x%016llx%016llx\n",nm,(unsigned long long)a.v[1],(unsigned long long)a.v[0]);}
static void q_to_mpfr(mpfr_t o,q a){uint64_t hi=a.v[1],lo=a.v[0]; int s=hi>>63; int e=(hi>>48)&0x7fff;
  uint64_t mh=hi&0xffffffffffffULL; mpfr_set_ui(o,mh,MPFR_RNDN); mpfr_mul_2ui(o,o,64,MPFR_RNDN);
  mpfr_t t; mpfr_init2(t,128); mpfr_set_ui(t,lo,MPFR_RNDN); mpfr_add(o,o,t,MPFR_RNDN); mpfr_clear(t);
  mpfr_t one; mpfr_init2(one,128); mpfr_set_ui(one,1,MPFR_RNDN); mpfr_mul_2ui(one,one,112,MPFR_RNDN);
  mpfr_add(o,o,one,MPFR_RNDN); mpfr_clear(one); mpfr_mul_2si(o,o,e-16383-112,MPFR_RNDN); if(s)mpfr_neg(o,o,MPFR_RNDN);}
static double ulp_err(q g,mpfr_t tr){
  mpfr_t gv,d,u; mpfr_inits2(300,gv,d,u,(mpfr_ptr)0); q_to_mpfr(gv,g);
  mpfr_sub(d,gv,tr,MPFR_RNDN); long e2; mpfr_get_d_2exp(&e2,gv,MPFR_RNDN);
  mpfr_set_ui(u,1,MPFR_RNDN); mpfr_mul_2si(u,u,(int)e2-1-112,MPFR_RNDN); mpfr_div(d,d,u,MPFR_RNDN);
  double r=mpfr_get_d(d,MPFR_RNDN); if(r<0)r=-r; mpfr_clears(gv,d,u,(mpfr_ptr)0); return r;}

static const q LOG2E={{0xe1777d0ffda0d23aull,0x3fff71547652b82full}},LN2HI={{0xf35793c600000000ull,0x3ffe62e42fefa39eull}},LN2LO={{0x81e6864ce5316c5bull,0x3fae673007e5ed5eull}},HALF={{0x0000000000000000ull,0x3ffe000000000000ull}},ONE={{0x0000000000000000ull,0x3fff000000000000ull}},SQRT2={{0xc908b2fb1366ea95ull,0x3fff6a09e667f3bcull}};
static const q EXC[25]={
  {{0x0000000000000000ull,0x3fff000000000000ull}}, {{0x0000000000000000ull,0x3fff000000000000ull}}, {{0x0000000000000000ull,0x3ffe000000000000ull}}, {{0x5555555555555555ull,0x3ffc555555555555ull}}, {{0x5555555555555555ull,0x3ffa555555555555ull}}, {{0x1111111111111111ull,0x3ff8111111111111ull}}, {{0x6c16c16c16c16c17ull,0x3ff56c16c16c16c1ull}}, {{0xa01a01a01a01a3e8ull,0x3ff2a01a01a01a01ull}}, {{0xa01a01a01a01a146ull,0x3fefa01a01a01a01ull}}, {{0x38faac1c88a5a526ull,0x3fec71de3a556c73ull}}, {{0xc72ef016d3d6e867ull,0x3fe927e4fb7789f5ull}}, {{0x38fe748363c46e8bull,0x3fe5ae64567f544eull}}, {{0x7b544dab18f475c5ull,0x3fe21eed8eff8d89ull}}, {{0x97c9f3aebabb2423ull,0x3fde6124613a86d0ull}}, {{0xd20b83c7f94d17d8ull,0x3fda93974a8c07c9ull}}, {{0xf5f4284f0d74f9e7ull,0x3fd6ae7f3e733b81ull}}, {{0xf417b4d27c5f92a9ull,0x3fd2ae7f3e733b81ull}}, {{0x6a419e674779c97cull,0x3fce952c77030a99ull}}, {{0x0466ff8c8b42b3dfull,0x3fca6827863b97b5ull}}, {{0x874b7a686d819241ull,0x3fc62f49b469f892ull}}, {{0xbb3b32a11bb5f139ull,0x3fc1e542ba427463ull}}, {{0xc93890ff9ab55cbbull,0x3fbd71b8db9f7f73ull}}, {{0x6efc0717eae785a1ull,0x3fb90ce38aab7bd7ull}}, {{0xcb3f4f7edfaa2666ull,0x3fb47693274bab2aull}}, {{0x61cb0e23655d47cbull,0x3faff3629154e0a7ull}}
};
static const q LOGC[23]={
  {{0x5555555555555555ull,0x3ffd555555555555ull}}, {{0x999999999999999aull,0x3ffc999999999999ull}}, {{0x2492492492492492ull,0x3ffc249249249249ull}}, {{0xc71c71c71c71c71cull,0x3ffbc71c71c71c71ull}}, {{0x5d1745d1745d1746ull,0x3ffb745d1745d174ull}}, {{0x3b13b13b13b13b14ull,0x3ffb3b13b13b13b1ull}}, {{0x1111111111111111ull,0x3ffb111111111111ull}}, {{0xe1e1e1e1e1e1e1e2ull,0x3ffae1e1e1e1e1e1ull}}, {{0x86bca1af286bca1bull,0x3ffaaf286bca1af2ull}}, {{0x8618618618618618ull,0x3ffa861861861861ull}}, {{0x42c8590b21642c86ull,0x3ffa642c8590b216ull}}, {{0xae147ae147ae147bull,0x3ffa47ae147ae147ull}}, {{0x84bda12f684bda13ull,0x3ffa2f684bda12f6ull}}, {{0x611a7b9611a7b961ull,0x3ffa1a7b9611a7b9ull}}, {{0x4210842108421084ull,0x3ffa084210842108ull}}, {{0x7c1f07c1f07c1f08ull,0x3ff9f07c1f07c1f0ull}}, {{0xd41d41d41d41d41dull,0x3ff9d41d41d41d41ull}}, {{0xf914c1bacf914c1cull,0x3ff9bacf914c1bacull}}, {{0xa41a41a41a41a41aull,0x3ff9a41a41a41a41ull}}, {{0x9c18f9c18f9c18faull,0x3ff98f9c18f9c18full}}, {{0x417d05f417d05f41ull,0x3ff97d05f417d05full}}, {{0x6c16c16c16c16c17ull,0x3ff96c16c16c16c1ull}}, {{0x72620ae4c415c988ull,0x3ff95c9882b93105ull}}
};
static q expq(q x){
  q t=qadd(qmul(x,LOG2E),HALF); int64_t k=f128M_to_i64(&t,softfloat_round_min,false); q qk=qi(k);
  q r=qsub(qsub(x,qmul(qk,LN2HI)),qmul(qk,LN2LO));
  q p=EXC[24]; for(int i=23;i>=0;i--) p=qadd(qmul(p,r),EXC[i]);
  return qmul(p,q2k((int)k));}
static q logq(q x){
  uint64_t hi=x.v[1]; int ef=((hi>>48)&0x7fff)-16383;
  q m; m.v[0]=x.v[0]; m.v[1]=(hi&0xffffffffffffULL)|(16383ULL<<48);
  if(f128M_le(&SQRT2,&m)){ m=qmul(m,HALF); ef++; }
  q f=qsub(m,ONE), s=qdiv(f,qadd(m,ONE)), z=qmul(s,s);
  q p2=LOGC[22]; for(int i=21;i>=0;i--) p2=qadd(qmul(p2,z),LOGC[i]);
  q r=qmul(qadd(z,z),p2), l1=qsub(f,qmul(s,qsub(f,r))), e=qi(ef);
  return qadd(qmul(e,LN2HI),qadd(l1,qmul(e,LN2LO)));}
static void sweep(const char*nm,q(*fn)(q),mpfr_t(*tru),int lo,int hi,int den){
  (void)tru; }
int main(void){
  softfloat_roundingMode=softfloat_round_near_even;
  mpfr_t tr; mpfr_init2(tr,300); double we=0,wl=0; q inv128=q2k(-7);
  for(int i=-2560;i<=2560;i++){ q x=qmul(qi(i),inv128); q g=expq(x);
    mpfr_set_si(tr,i,MPFR_RNDN); mpfr_div_ui(tr,tr,128,MPFR_RNDN); mpfr_exp(tr,tr,MPFR_RNDN);
    double e=ulp_err(g,tr); if(e>we)we=e; }
  for(int i=1;i<=25600;i++){ q x=qmul(qi(i),inv128); q g=logq(x);
    mpfr_set_si(tr,i,MPFR_RNDN); mpfr_div_ui(tr,tr,128,MPFR_RNDN); mpfr_log(tr,tr,MPFR_RNDN); if(i==128)continue;
    double e=ulp_err(g,tr); if(e>wl)wl=e; }
  printf("# rq exp max %.3f ULP [-20,20];  log max %.3f ULP (0,200]\n",we,wl);
  qprint("exp(1)  = ",expq(qi(1)));
  qprint("log(2)  = ",logq(qi(2)));
  qprint("log(10) = ",logq(qi(10)));
  mpfr_clear(tr); return 0;}
