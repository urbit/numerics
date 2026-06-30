// Standalone validation of the native arithmetic in twoc.c.
// Copies the _TWOC_NATIVE macro + mask helpers VERBATIM from the jet, and
// checks them against an independent two's-complement reference (u128/s128).
#include <stdio.h>
#include <stdlib.h>

typedef unsigned long long c3_d;
typedef unsigned char      c3_t;
typedef unsigned __int128  _twoc_u128;
typedef          __int128  _twoc_s128;
#define c3y 0
#define c3n 1

// ---- VERBATIM from twoc.c ----
#define _TWOC_NATIVE(SUF, UT)                                                  \
  static inline UT _tn_neg##SUF(UT a, UT msk) {                                \
    return (UT)(((~a) + 1) & msk);                                             \
  }                                                                            \
  static inline UT _tn_sign##SUF(UT a, c3_d sb) {                              \
    return (UT)((a >> sb) & 1);                                               \
  }                                                                            \
  static inline UT _tn_abs##SUF(UT a, UT msk, c3_d sb) {                       \
    return _tn_sign##SUF(a, sb) ? _tn_neg##SUF(a, msk) : (UT)(a & msk);        \
  }                                                                            \
  static inline UT _tn_add##SUF(UT a, UT b, UT msk) {                          \
    return (UT)((a + b) & msk);                                                \
  }                                                                            \
  static inline UT _tn_sub##SUF(UT a, UT b, UT msk) {                          \
    return (UT)((a + _tn_neg##SUF(b, msk)) & msk);                             \
  }                                                                            \
  static inline UT _tn_mul##SUF(UT a, UT b, UT msk) {                          \
    return (UT)(((a & msk) * (b & msk)) & msk);                                \
  }                                                                            \
  static inline UT _tn_div##SUF(UT a, UT b, UT msk, c3_d sb) {                 \
    UT qa = _tn_abs##SUF(a, msk, sb), qb = _tn_abs##SUF(b, msk, sb);           \
    UT q  = (UT)(qa / qb);                                                     \
    return ( _tn_sign##SUF(a, sb) != _tn_sign##SUF(b, sb) )                    \
           ? _tn_neg##SUF(q, msk) : (UT)(q & msk);                            \
  }                                                                            \
  static inline UT _tn_rem##SUF(UT a, UT b, UT msk, c3_d sb) {                 \
    UT r = (UT)(_tn_abs##SUF(a, msk, sb) % _tn_abs##SUF(b, msk, sb));          \
    return _tn_sign##SUF(a, sb) ? _tn_neg##SUF(r, msk) : (UT)(r & msk);        \
  }                                                                            \
  static inline UT _tn_pow##SUF(UT a, c3_d n, UT msk) {                        \
    UT base = (UT)(a & msk), acc = (UT)(1 & msk);                              \
    while ( n ) {                                                              \
      if ( n & 1 ) acc = (UT)((acc * base) & msk);                            \
      base = (UT)((base * base) & msk);                                        \
      n >>= 1;                                                                 \
    }                                                                          \
    return acc;                                                               \
  }                                                                            \
  static inline c3_t _tn_gth##SUF(UT a, UT b, c3_d sb) {                       \
    UT sa = _tn_sign##SUF(a, sb), sbb = _tn_sign##SUF(b, sb);                  \
    if ( sa != sbb ) return ( 0 == sa ) ? c3y : c3n;                          \
    return ( a > b ) ? c3y : c3n;                                             \
  }

_TWOC_NATIVE(64,  c3_d)
_TWOC_NATIVE(128, _twoc_u128)

static inline c3_d _tn_msk64(c3_d wid) {
  return ( wid >= 64 ) ? (c3_d)~0ull : (((c3_d)1 << wid) - 1);
}
static inline _twoc_u128 _tn_msk128(c3_d wid) {
  return ( wid >= 128 ) ? (~(_twoc_u128)0) : ((((_twoc_u128)1) << wid) - 1);
}
static inline c3_t _twoc_not(c3_t b) { return ( c3y == b ) ? c3n : c3y; }
// ---- end verbatim ----

typedef unsigned __int128 u128;
typedef          __int128 s128;

// reference (wid <= 100 so modv fits u128 and signed fits s128)
static u128 R_mask(int w){ return (((u128)1) << w) - 1; }
static s128 R_sgn(u128 p, int w){
  u128 m = R_mask(w);
  p &= m;
  if ( (p >> (w-1)) & 1 ) return (s128)p - (s128)(((u128)1) << w);
  return (s128)p;
}
static u128 R_pat(s128 v, int w){ return ((u128)v) & R_mask(w); }

static long long FAILS = 0, TOTS = 0;
static void chk(const char* op, int w, u128 a, u128 b, u128 got, u128 ref){
  TOTS++;
  if ( got != ref ) {
    FAILS++;
    if ( FAILS <= 40 ) {
      fprintf(stderr, "FAIL %s w=%d a=%016llx%016llx b=%016llx%016llx got=%016llx%016llx ref=%016llx%016llx\n",
        op, w, (c3_d)(a>>64),(c3_d)a, (c3_d)(b>>64),(c3_d)b,
        (c3_d)(got>>64),(c3_d)got, (c3_d)(ref>>64),(c3_d)ref);
    }
  }
}

// dispatch wrappers
static u128 J_add(u128 a,u128 b,int w){ return w<=64? _tn_add64(a,b,_tn_msk64(w)) : _tn_add128(a,b,_tn_msk128(w)); }
static u128 J_sub(u128 a,u128 b,int w){ return w<=64? _tn_sub64(a,b,_tn_msk64(w)) : _tn_sub128(a,b,_tn_msk128(w)); }
static u128 J_mul(u128 a,u128 b,int w){ return w<=64? _tn_mul64(a,b,_tn_msk64(w)) : _tn_mul128(a,b,_tn_msk128(w)); }
static u128 J_div(u128 a,u128 b,int w){ return w<=64? _tn_div64(a,b,_tn_msk64(w),w-1) : _tn_div128(a,b,_tn_msk128(w),w-1); }
static u128 J_rem(u128 a,u128 b,int w){ return w<=64? _tn_rem64(a,b,_tn_msk64(w),w-1) : _tn_rem128(a,b,_tn_msk128(w),w-1); }
static u128 J_neg(u128 a,int w){ return w<=64? _tn_neg64(a,_tn_msk64(w)) : _tn_neg128(a,_tn_msk128(w)); }
static u128 J_abs(u128 a,int w){ return w<=64? _tn_abs64(a,_tn_msk64(w),w-1) : _tn_abs128(a,_tn_msk128(w),w-1); }
static c3_t J_gth(u128 a,u128 b,int w){ return w<=64? _tn_gth64(a,b,w-1) : _tn_gth128(a,b,w-1); }
static u128 J_pow(u128 a,c3_d n,int w){ return w<=64? _tn_pow64(a,n,_tn_msk64(w)) : _tn_pow128(a,n,_tn_msk128(w)); }

int main(void){
  int widths[] = {8,16,17,32,64,100};
  for ( int wi=0; wi<6; wi++ ){
    int w = widths[wi];
    u128 m = R_mask(w);
    u128 half = ((u128)1) << (w-1);
    u128 vals[] = { 0,1,2,3, half-1, half, half+1, m-1, m, (u128)5<<(w/2), m/3, (m/7)|half };
    int NV = sizeof(vals)/sizeof(vals[0]);
    for ( int i=0;i<NV;i++ ) for ( int j=0;j<NV;j++ ){
      u128 a=vals[i]&m, b=vals[j]&m;
      // add/sub/mul
      chk("add",w,a,b, J_add(a,b,w), (a+b)&m);
      chk("sub",w,a,b, J_sub(a,b,w), R_pat(R_sgn(a,w)-R_sgn(b,w),w));
      chk("mul",w,a,b, J_mul(a,b,w), (a*b)&m);
      // gth as 0/1 mapped: J returns c3y(0)/c3n(1); ref bool
      c3_t jg = J_gth(a,b,w);
      u128 rg = (R_sgn(a,w) > R_sgn(b,w)) ? c3y : c3n;
      chk("gth",w,a,b, jg, rg);
      c3_t jl = J_gth(b,a,w);
      chk("lth",w,a,b, jl, (R_sgn(a,w) < R_sgn(b,w))?c3y:c3n);
      chk("lte",w,a,b, _twoc_not(J_gth(a,b,w)), (R_sgn(a,w) <= R_sgn(b,w))?c3y:c3n);
      chk("gte",w,a,b, _twoc_not(J_gth(b,a,w)), (R_sgn(a,w) >= R_sgn(b,w))?c3y:c3n);
      // div/rem (skip zero divisor)
      if ( (b&m) != 0 ){
        s128 sa=R_sgn(a,w), sb=R_sgn(b,w);
        s128 q = sa / sb;            // trunc toward zero
        s128 r = sa - q*sb;
        chk("div",w,a,b, J_div(a,b,w), R_pat(q,w));
        chk("rem",w,a,b, J_rem(a,b,w), R_pat(r,w));
      }
    }
    // unary + pow
    for ( int i=0;i<NV;i++ ){
      u128 a=vals[i]&m;
      chk("neg",w,a,0, J_neg(a,w), R_pat(-R_sgn(a,w),w));
      chk("abs",w,a,0, J_abs(a,w), (R_sgn(a,w)<0)? R_pat(-R_sgn(a,w),w) : (a&m));
      c3_d ns[]={0,1,2,3,5,8,13};
      for ( int k=0;k<7;k++ ){
        c3_d n=ns[k];
        // independent reference pow: naive repeated multiply, a^n mod 2^w
        u128 acc=1&m, base=a&m;
        for ( c3_d t=0; t<n; t++ ) acc=(acc*base)&m;
        chk("pow",w,a,n, J_pow(a,n,w), acc);
      }
    }
  }
  printf("checks=%lld fails=%lld\n", TOTS, FAILS);
  return FAILS ? 1 : 0;
}
