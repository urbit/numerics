// rq_check.c -- SoftFloat-f128 reference for @rq transcendentals (PR #18).
// SoftFloat = the jet arithmetic (this code IS the jet body, modulo the u3 ABI);
// MPFR = correctly-rounded truth.  Coeffs from mpmath @113-bit.
// Build: re-archive libsoftfloat.a with `libtool -static` (zig .a not 8-byte
// aligned for Apple ld), then cc -I<sf-include> -I/opt/homebrew/include
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
static q qsqt(q a){q r; f128M_sqrt(&a,&r); return r;}
static q qi(int64_t n){q r; i64_to_f128M(n,&r); return r;}
static q q2k(int k){q r; r.v[0]=0; r.v[1]=((uint64_t)(k+16383))<<48; return r;}
static q negq(q a){a.v[1]^=0x8000000000000000ULL; return a;}
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
static const q LOG2E={{0xe1777d0ffda0d23aull,0x3fff71547652b82full}},LN2HI={{0xf35793c800000000ull,0x3ffe62e42fefa39eull}},LN2LO={{0xfc32f366359d274aull,0xbfad319ff0342542ull}},HALF={{0x0000000000000000ull,0x3ffe000000000000ull}},ONE={{0x0000000000000000ull,0x3fff000000000000ull}},TWO={{0x0000000000000000ull,0x4000000000000000ull}},SQRT2={{0xc908b2fb1366ea95ull,0x3fff6a09e667f3bcull}};
static const q INVPIO2={{0x2a53f84eafa3ea6aull,0x3ffe45f306dc9c88ull}},PIO2_1={{0x8460000000000000ull,0x3fff921fb54442d1ull}},PIO2_1T={{0x07344a409382229aull,0x3fc2313198a2e037ull}};
static const q THIRD={{0x5555555555555555ull,0x3ffd555555555555ull}};  // 1/3, matches ++cbt
// exp: fdlibm rational reconstruction; EXC = even minimax P(t), t=r^2 (deg-10)
static const q EXC[11]={
  {{0x5555555555555555ull,0x3ffc555555555555ull}},
  {{0x6c16c16c16c09e83ull,0xbff66c16c16c16c1ull}},
  {{0x6abc0115453d96ddull,0x3ff11566abc01156ull}},
  {{0xaac663e4a6d65ccaull,0xbfebbbd779334ef0ull}},
  {{0xda06115986f507fbull,0x3fe666a8f2bf70ebull}},
  {{0x43eb0e288c2e45a8ull,0xbfe122805d644267ull}},
  {{0x12be0476b628552full,0x3fdbd6db2c4e0507ull}},
  {{0xeb838f5da821635aull,0xbfd67da4e1efb419ull}},
  {{0xdc61daecbfc0d781ull,0x3fd1355867f7df64ull}},
  {{0x54bb7852bc52bd9aull,0xbfcbf56e4264f8adull}},
  {{0x822162270789ca71ull,0x3fc68fc13579bfe0ull}}
};
static const q LOGC[23]={
  {{0x5555555555555555ull,0x3ffd555555555555ull}}, {{0x999999999999999aull,0x3ffc999999999999ull}}, {{0x2492492492492492ull,0x3ffc249249249249ull}}, {{0xc71c71c71c71c71cull,0x3ffbc71c71c71c71ull}}, {{0x5d1745d1745d1746ull,0x3ffb745d1745d174ull}}, {{0x3b13b13b13b13b14ull,0x3ffb3b13b13b13b1ull}}, {{0x1111111111111111ull,0x3ffb111111111111ull}}, {{0xe1e1e1e1e1e1e1e2ull,0x3ffae1e1e1e1e1e1ull}}, {{0x86bca1af286bca1bull,0x3ffaaf286bca1af2ull}}, {{0x8618618618618618ull,0x3ffa861861861861ull}}, {{0x42c8590b21642c86ull,0x3ffa642c8590b216ull}}, {{0xae147ae147ae147bull,0x3ffa47ae147ae147ull}}, {{0x84bda12f684bda13ull,0x3ffa2f684bda12f6ull}}, {{0x611a7b9611a7b961ull,0x3ffa1a7b9611a7b9ull}}, {{0x4210842108421084ull,0x3ffa084210842108ull}}, {{0x7c1f07c1f07c1f08ull,0x3ff9f07c1f07c1f0ull}}, {{0xd41d41d41d41d41dull,0x3ff9d41d41d41d41ull}}, {{0xf914c1bacf914c1cull,0x3ff9bacf914c1bacull}}, {{0xa41a41a41a41a41aull,0x3ff9a41a41a41a41ull}}, {{0x9c18f9c18f9c18faull,0x3ff98f9c18f9c18full}}, {{0x417d05f417d05f41ull,0x3ff97d05f417d05full}}, {{0x6c16c16c16c16c17ull,0x3ff96c16c16c16c1ull}}, {{0x72620ae4c415c988ull,0x3ff95c9882b93105ull}}
};
static const q SINC[16]={
  {{0x5555555555555555ull,0xbffc555555555555ull}}, {{0x1111111111111111ull,0x3ff8111111111111ull}}, {{0xa01a01a01a01a01aull,0xbff2a01a01a01a01ull}}, {{0x38faac1c88e50017ull,0x3fec71de3a556c73ull}}, {{0x38fe747e4b837dc7ull,0xbfe5ae64567f544eull}}, {{0x97ca38331d23af68ull,0x3fde6124613a86d0ull}}, {{0xf11d8656b0ee8cb0ull,0xbfd6ae7f3e733b81ull}}, {{0xa6b2605197771b00ull,0x3fce952c77030ad4ull}}, {{0x724ca1ec3b7b9675ull,0xbfc62f49b4681415ull}}, {{0x18bef146fcee6e45ull,0x3fbd71b8ef6dcf57ull}}, {{0x9d97b8704dd7f628ull,0xbfb4761b41316381ull}}, {{0x8d4e44a419776f11ull,0x3fab3f3ccdd165faull}}, {{0x320a9a18f15d4277ull,0xbfa1d1ab1c2dcceaull}}, {{0xd7abe30e7766f129ull,0x3f98259f98b4358aull}}, {{0xc42e1ee46fa6bfc4ull,0xbf8e434d2e783f5bull}}, {{0x1b5382cdffa97422ull,0x3f843981254dd0d5ull}}
};
static const q COSC[16]={
  {{0x5555555555555555ull,0x3ffa555555555555ull}}, {{0x6c16c16c16c16c17ull,0xbff56c16c16c16c1ull}}, {{0xa01a01a01a01a01aull,0x3fefa01a01a01a01ull}}, {{0xc72ef016d3ea6679ull,0xbfe927e4fb7789f5ull}}, {{0x7b544da987acfe85ull,0x3fe21eed8eff8d89ull}}, {{0xd20badf145dfa3e5ull,0xbfda93974a8c07c9ull}}, {{0xf11d8656b0ee8cb0ull,0x3fd2ae7f3e733b81ull}}, {{0x77bb004886a2c2abull,0xbfca6827863b97d9ull}}, {{0x507a9cad2bf8f0bbull,0x3fc1e542ba402022ull}}, {{0x29450c90b7f338ecull,0xbfb90ce396db7f85ull}}, {{0x7cca4b4067ca9d8aull,0x3faff2cf01972f57ull}}, {{0x9a38f2050ba6b015ull,0xbfa688e85fc6a4e5ull}}, {{0xd373c5c51c354a8dull,0x3f9d0a18a2635085ull}}, {{0xe60caded4c2989c5ull,0xbf933932c5047d60ull}}, {{0xc42e1ee46fa6bfc4ull,0x3f89434d2e783f5bull}}, {{0xa13f8a2b4af9d6b7ull,0xbf7f2710231c0fd7ull}}
};
static q expq(q x){
  q t0=qadd(qmul(x,LOG2E),HALF); int64_t k=f128M_to_i64(&t0,softfloat_round_min,false); q qk=qi(k);
  q hi=qsub(x,qmul(qk,LN2HI)), lo=qmul(qk,LN2LO), r=qsub(hi,lo);
  q t=qmul(r,r);
  q c=EXC[10]; for(int i=9;i>=0;i--) c=qadd(qmul(c,t),EXC[i]);   // P(t)
  c=qsub(r,qmul(t,c));                                           // c = r - t*P
  q y=qsub(ONE, qsub(qsub(lo, qdiv(qmul(r,c),qsub(TWO,c))), hi));// 1-((lo-r*c/(2-c))-hi)
  return qmul(y,q2k((int)k));}
static q logq(q x){
  uint64_t hi=x.v[1]; int ef=((hi>>48)&0x7fff)-16383;
  q m; m.v[0]=x.v[0]; m.v[1]=(hi&0xffffffffffffULL)|(16383ULL<<48);
  if(f128M_le(&SQRT2,&m)){ m=qmul(m,HALF); ef++; }
  q f=qsub(m,ONE), s=qdiv(f,qadd(m,ONE)), z=qmul(s,s);
  q p2=LOGC[22]; for(int i=21;i>=0;i--) p2=qadd(qmul(p2,z),LOGC[i]);
  q r=qmul(qadd(z,z),p2), l1=qsub(f,qmul(s,qsub(f,r))), e=qi(ef);
  return qadd(qmul(e,LN2HI),qadd(l1,qmul(e,LN2LO)));}
static q ksinq(q x,q y){
  q z=qmul(x,x);
  q r=SINC[15]; for(int i=14;i>=1;i--) r=qadd(qmul(r,z),SINC[i]);
  q v=qmul(z,x), aa=qsub(qmul(HALF,y),qmul(v,r)), bb=qsub(qmul(z,aa),y), cc=qsub(bb,qmul(v,SINC[0]));
  return qsub(x,cc);}
static q kcosq(q x,q y){
  q z=qmul(x,x);
  q rc=COSC[15]; for(int i=14;i>=0;i--) rc=qadd(qmul(rc,z),COSC[i]);
  q hz=qmul(HALF,z), w2=qsub(ONE,hz);
  q aa=qsub(qsub(ONE,w2),hz), bb=qsub(qmul(qmul(z,z),rc),qmul(x,y));
  return qadd(w2,qadd(aa,bb));}
static void reduce(q ax,int64_t*qq,q*rhi,q*rlo){
  q t=qadd(qmul(ax,INVPIO2),HALF); int64_t k=f128M_to_i64(&t,softfloat_round_min,false);
  q qk=qi(k), hp=qsub(ax,qmul(qk,PIO2_1)), w=qmul(qk,PIO2_1T);
  *qq=k; *rhi=qsub(hp,w); *rlo=qsub(qsub(hp,*rhi),w);}
static q sinq(q x){
  int neg=x.v[1]>>63; q ax=x; ax.v[1]&=0x7fffffffffffffffULL;
  int64_t k; q rhi,rlo; reduce(ax,&k,&rhi,&rlo); int m=((int)k)&3;
  q ks=ksinq(rhi,rlo), kc=kcosq(rhi,rlo);
  q v = m==0?ks : m==1?kc : m==2?negq(ks):negq(kc);
  return neg?negq(v):v;}
static q cosq(q x){
  q ax=x; ax.v[1]&=0x7fffffffffffffffULL;
  int64_t k; q rhi,rlo; reduce(ax,&k,&rhi,&rlo); int m=((int)k)&3;
  q ks=ksinq(rhi,rlo), kc=kcosq(rhi,rlo);
  return m==0?kc : m==1?negq(ks) : m==2?negq(kc):ks;}
static q cbtq(q x){  // ++cbt: exp(log|x|/3), sign-restored (perfect cubes -> exact)
  q ax=x; ax.v[1]&=0x7fffffffffffffffULL;
  q r=expq(qmul(logq(ax),THIRD));
  return (x.v[1]>>63) ? negq(r) : r;}
int main(void){
  softfloat_roundingMode=softfloat_round_near_even;
  mpfr_t tr; mpfr_init2(tr,300); double we=0,wl=0,ws=0,wc=0,wq=0; q inv128=q2k(-7);
  for(int i=-2560;i<=2560;i++){ q x=qmul(qi(i),inv128); q g=expq(x);
    mpfr_set_si(tr,i,MPFR_RNDN); mpfr_div_ui(tr,tr,128,MPFR_RNDN); mpfr_exp(tr,tr,MPFR_RNDN);
    double e=ulp_err(g,tr); if(e>we)we=e; }
  for(int i=1;i<=25600;i++){ q x=qmul(qi(i),inv128); q g=logq(x);
    mpfr_set_si(tr,i,MPFR_RNDN); mpfr_div_ui(tr,tr,128,MPFR_RNDN); mpfr_log(tr,tr,MPFR_RNDN); if(i==128)continue;
    double e=ulp_err(g,tr); if(e>wl)wl=e; }
  for(int i=-12800;i<=12800;i++){ q x=qmul(qi(i),inv128);
    mpfr_set_si(tr,i,MPFR_RNDN); mpfr_div_ui(tr,tr,128,MPFR_RNDN);
    mpfr_t ts,tc; mpfr_init2(ts,300); mpfr_init2(tc,300); mpfr_sin(ts,tr,MPFR_RNDN); mpfr_cos(tc,tr,MPFR_RNDN);
    double es=mpfr_get_d(ts,MPFR_RNDN), ec=mpfr_get_d(tc,MPFR_RNDN);
    if(es>1e-9||es<-1e-9){double e=ulp_err(sinq(x),ts); if(e>ws)ws=e;}
    if(ec>1e-9||ec<-1e-9){double e=ulp_err(cosq(x),tc); if(e>wc)wc=e;}
    mpfr_clear(ts); mpfr_clear(tc); }
  for(int i=1;i<=25600;i++){ q x=qmul(qi(i),inv128); q g=qsqt(x);
    mpfr_set_si(tr,i,MPFR_RNDN); mpfr_div_ui(tr,tr,128,MPFR_RNDN); mpfr_sqrt(tr,tr,MPFR_RNDN);
    double e=ulp_err(g,tr); if(e>wq)wq=e; }
  double wb=0;                                            // cbt = exp(log|x|/3)
  for(int i=-25600;i<=25600;i++){ if(i==0)continue; q x=qmul(qi(i),inv128); q g=cbtq(x);
    mpfr_set_si(tr,i,MPFR_RNDN); mpfr_div_ui(tr,tr,128,MPFR_RNDN); mpfr_cbrt(tr,tr,MPFR_RNDN);
    double e=ulp_err(g,tr); if(e>wb)wb=e; }
  printf("# rq: exp %.3f  log %.3f  sin %.3f  cos %.3f  sqrt %.3f  cbt %.3f ULP\n",we,wl,ws,wc,wq,wb);
  qprint("sin(1)   = ",sinq(qi(1)));
  qprint("cos(1)   = ",cosq(qi(1)));
  qprint("sqrt(2)  = ",qsqt(qi(2)));
  qprint("cbt(-27) = ",cbtq(qi(-27)));               // perfect cube -> exact -3 = 0xc000.8000...
  mpfr_clear(tr); return 0;}
