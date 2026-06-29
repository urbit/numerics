// Validate the GMP (arbitrary-precision) path of twoc.c.
// The _tg_* helpers + op blocks are copied VERBATIM from the jet; we run them
// at widths 8..100 and compare to the independent hardware u128/s128 reference
// (the same one that validated the native path), then hand-check a few >128
// cases that no hardware type can hold.
#include <stdio.h>
#include <gmp.h>

typedef unsigned long long c3_d;
typedef unsigned char      c3_t;
typedef unsigned __int128  u128;
typedef          __int128  s128;
#define c3y 0
#define c3n 1

// ---- VERBATIM helpers from twoc.c ----
static inline void _tg_mask(mpz_t r, const mpz_t a, c3_d wid){ mpz_fdiv_r_2exp(r,a,wid); }
static inline c3_t _tg_sign(const mpz_t m, c3_d wid){ return mpz_tstbit(m,wid-1)?c3y:c3n; }
static inline void _tg_signed(mpz_t s, const mpz_t m, c3_d wid){
  if ( _tg_sign(m,wid)==c3y ){ mpz_t p; mpz_init(p); mpz_ui_pow_ui(p,2,wid);
    mpz_sub(s,m,p); mpz_clear(p);
  } else mpz_set(s,m);
}
// ---- op blocks mirroring the jet GMPBLOCKs (operate on in-range inputs) ----
static void G_add(mpz_t r,mpz_t a,mpz_t b,c3_d w){ mpz_add(r,a,b); _tg_mask(r,r,w); }
static void G_sub(mpz_t r,mpz_t a,mpz_t b,c3_d w){ mpz_sub(r,a,b); _tg_mask(r,r,w); }
static void G_mul(mpz_t r,mpz_t a,mpz_t b,c3_d w){ mpz_mul(r,a,b); _tg_mask(r,r,w); }
static void G_neg(mpz_t r,mpz_t a,c3_d w){ mpz_neg(r,a); _tg_mask(r,r,w); }
static void G_abs(mpz_t r,mpz_t a,c3_d w){ mpz_t m; mpz_init(m); _tg_mask(m,a,w);
  if(_tg_sign(m,w)==c3y){ mpz_neg(m,m); _tg_mask(m,m,w);} mpz_set(r,m); mpz_clear(m); }
static void G_div(mpz_t r,mpz_t a,mpz_t b,c3_d w){ mpz_t am,bm,as,bs;
  mpz_inits(am,bm,as,bs,NULL); _tg_mask(am,a,w); _tg_mask(bm,b,w);
  _tg_signed(as,am,w); _tg_signed(bs,bm,w); mpz_tdiv_q(as,as,bs); _tg_mask(as,as,w);
  mpz_set(r,as); mpz_clears(am,bm,as,bs,NULL); }
static void G_rem(mpz_t r,mpz_t a,mpz_t b,c3_d w){ mpz_t am,bm,as,bs;
  mpz_inits(am,bm,as,bs,NULL); _tg_mask(am,a,w); _tg_mask(bm,b,w);
  _tg_signed(as,am,w); _tg_signed(bs,bm,w); mpz_tdiv_r(as,as,bs); _tg_mask(as,as,w);
  mpz_set(r,as); mpz_clears(am,bm,as,bs,NULL); }
static void G_pow(mpz_t r,mpz_t a,c3_d n,c3_d w){ mpz_t am,nm,mod;
  mpz_inits(am,nm,mod,NULL); _tg_mask(am,a,w); mpz_set_ui(nm,n);
  mpz_ui_pow_ui(mod,2,w); mpz_powm(am,am,nm,mod); mpz_set(r,am);
  mpz_clears(am,nm,mod,NULL); }
static c3_t G_gth(mpz_t a,mpz_t b,c3_d w){ mpz_t as,bs; mpz_inits(as,bs,NULL);
  _tg_signed(as,a,w); _tg_signed(bs,b,w); int c=mpz_cmp(as,bs);
  mpz_clears(as,bs,NULL); return c>0?c3y:c3n; }

// ---- independent hardware reference (wid <= 100) ----
static u128 R_mask(int w){ return (((u128)1)<<w)-1; }
static s128 R_sgn(u128 p,int w){ p&=R_mask(w);
  if((p>>(w-1))&1) return (s128)p-(s128)(((u128)1)<<w); return (s128)p; }
static u128 R_pat(s128 v,int w){ return ((u128)v)&R_mask(w); }

static void u128_to_mpz(mpz_t r,u128 v){ mpz_set_ui(r,(c3_d)(v>>64));
  mpz_mul_2exp(r,r,64); mpz_t lo; mpz_init_set_ui(lo,(c3_d)v); mpz_add(r,r,lo); mpz_clear(lo); }
static u128 mpz_to_u128(const mpz_t v){ mpz_t t,lo; mpz_inits(t,lo,NULL);
  mpz_fdiv_r_2exp(lo,v,64); mpz_fdiv_q_2exp(t,v,64);
  u128 r=((u128)mpz_get_ui(t)<<64)|(u128)mpz_get_ui(lo); mpz_clears(t,lo,NULL); return r; }

static long long FAILS=0, TOTS=0;
static void chku(const char*op,int w,u128 got,u128 ref){ TOTS++;
  if(got!=ref){ FAILS++; if(FAILS<=40)
    fprintf(stderr,"FAIL %s w=%d got=%llx%016llx ref=%llx%016llx\n",op,w,
      (c3_d)(got>>64),(c3_d)got,(c3_d)(ref>>64),(c3_d)ref); } }

int main(void){
  int widths[]={8,16,17,32,64,100};
  mpz_t a,b,r; mpz_inits(a,b,r,NULL);
  for(int wi=0;wi<6;wi++){ int w=widths[wi]; u128 m=R_mask(w), half=((u128)1)<<(w-1);
    u128 vals[]={0,1,2,3,half-1,half,half+1,m-1,m,(u128)5<<(w/2),m/3,(m/7)|half};
    int NV=12;
    for(int i=0;i<NV;i++)for(int j=0;j<NV;j++){
      u128 ua=vals[i]&m, ub=vals[j]&m; u128_to_mpz(a,ua); u128_to_mpz(b,ub);
      G_add(r,a,b,w); chku("add",w,mpz_to_u128(r),(ua+ub)&m);
      G_sub(r,a,b,w); chku("sub",w,mpz_to_u128(r),R_pat(R_sgn(ua,w)-R_sgn(ub,w),w));
      G_mul(r,a,b,w); chku("mul",w,mpz_to_u128(r),(ua*ub)&m);
      chku("gth",w,G_gth(a,b,w),(R_sgn(ua,w)>R_sgn(ub,w))?c3y:c3n);
      if((ub&m)!=0){ s128 sa=R_sgn(ua,w),sb=R_sgn(ub,w); s128 q=sa/sb, rr=sa-q*sb;
        G_div(r,a,b,w); chku("div",w,mpz_to_u128(r),R_pat(q,w));
        G_rem(r,a,b,w); chku("rem",w,mpz_to_u128(r),R_pat(rr,w)); }
    }
    for(int i=0;i<NV;i++){ u128 ua=vals[i]&m; u128_to_mpz(a,ua);
      G_neg(r,a,w); chku("neg",w,mpz_to_u128(r),R_pat(-R_sgn(ua,w),w));
      G_abs(r,a,w); chku("abs",w,mpz_to_u128(r),(R_sgn(ua,w)<0)?R_pat(-R_sgn(ua,w),w):(ua&m));
      c3_d ns[]={0,1,2,3,5,8,13};
      for(int k=0;k<7;k++){ c3_d n=ns[k]; u128 acc=1&m,base=ua&m;
        for(c3_d t=0;t<n;t++) acc=(acc*base)&m;
        G_pow(r,a,n,w); chku("pow",w,mpz_to_u128(r),acc); }
    }
  }
  // hand-checked wide cases (w=200): -1 + -1 = -2; -1 * -1 = 1; min/-1 = min; max+1 = min
  { c3_d w=200; mpz_t one,negone,res,exp; mpz_inits(one,negone,res,exp,NULL);
    mpz_set_ui(one,1); _tg_mask(negone,/*=-1 as 2^w-1*/ negone,w);
    mpz_ui_pow_ui(negone,2,w); mpz_sub_ui(negone,negone,1);   // 2^200 - 1 == -1
    G_add(res,negone,negone,w); mpz_ui_pow_ui(exp,2,w); mpz_sub_ui(exp,exp,2);
    TOTS++; if(mpz_cmp(res,exp)!=0){FAILS++; fprintf(stderr,"FAIL wide add -1+-1\n");}
    G_mul(res,negone,negone,w); TOTS++; if(mpz_cmp_ui(res,1)!=0){FAILS++; fprintf(stderr,"FAIL wide mul -1*-1\n");}
    mpz_t mn; mpz_init(mn); mpz_ui_pow_ui(mn,2,w-1);          // min = 2^199
    G_div(res,mn,negone,w); TOTS++; if(mpz_cmp(res,mn)!=0){FAILS++; fprintf(stderr,"FAIL wide div min/-1\n");}
    mpz_t mx; mpz_init(mx); mpz_ui_pow_ui(mx,2,w-1); mpz_sub_ui(mx,mx,1); // max=2^199-1
    G_add(res,mx,one,w); TOTS++; if(mpz_cmp(res,mn)!=0){FAILS++; fprintf(stderr,"FAIL wide max+1\n");}
    mpz_clears(one,negone,res,exp,mn,mx,NULL);
  }
  mpz_clears(a,b,r,NULL);
  printf("checks=%lld fails=%lld\n",TOTS,FAILS);
  return FAILS?1:0;
}
