/// @file
///
/// Jets for the numerics `math.hoon` library (userspace, registered under the
/// `non` chapter alongside `lagoon`).  Each body runs the IDENTICAL algorithm
/// as its Hoon arm (same reduction, coefficients, Horner order) in Berkeley
/// SoftFloat, so jet output is BIT-EXACT to the pure-Hoon reference -- not
/// merely faithful.  Transcendentals force round-nearest-even internally (they
/// take no rounding-mode axis), matching the Hoon.
///
/// Marshalling uses chub (64-bit) reads/writes so it is word-size-agnostic
/// across the 32-bit and 64-bit runtimes (the divergence that broke @rq sub).
///
/// MASTER COPY lives in urbit/numerics libmath/vere/; applied by hand to the
/// vere runtime (pkg/noun/jets/i/math.c) -- see libmath/vere/README.md.  The
/// `_rd_*` algorithm cores below are shared with the bit-exact harness
/// (libmath/vere/test/) via -DMATH_JET_HARNESS, so the jet and its test can
/// never drift.

#ifndef MATH_JET_HARNESS
#include "jets/q.h"
#include "jets/w.h"
#include "noun.h"
#endif
#include "softfloat.h"

#ifdef MATH_JET_HARNESS
#include <stdint.h>
typedef uint64_t c3_d;
typedef int64_t  c3_ds;
#endif

  union doub {
    float64_t d;
    c3_d c;
  };

  static const c3_d _RD_QNAN = 0x7ff8000000000000ULL;
  static const c3_d _RD_PINF = 0x7ff0000000000000ULL;
  static const c3_d _RD_NINF = 0xfff0000000000000ULL;

  static inline union doub _rd_bits(c3_d b) { union doub u; u.c = b; return u; }

/* @rd exp -- math.hoon ++rd ++exp
**   x = k*ln2 + r (Cody-Waite), exp(x) = 2^k * P(r), P a degree-11 minimax
**   polynomial; +scale2 is a correctly-rounded ldexp.
*/
  //  pow2(j) = (j+1023)<<52 as f64 bits = 2^j  (normal range j in [-1022,1023])
  static inline float64_t
  _rd_pow2(c3_ds j)
  {
    union doub u;
    u.c = ((c3_d)(j + 1023)) << 52;
    return u.d;
  }

  //  scale2: ldexp with overflow/subnormal tails (math.hoon:1519)
  static inline float64_t
  _rd_scale2(float64_t p, c3_ds k)
  {
    if ( (k - 1024) >= 0 ) {                         // k>=1024
      return f64_mul(f64_mul(p, _rd_pow2(1023)), _rd_pow2(k - 1023));
    }
    if ( !((k + 1022) >= 0) ) {                      // k<-1022
      return f64_mul(f64_mul(p, _rd_pow2(k + 54)), _rd_pow2(-54));
    }
    return f64_mul(p, _rd_pow2(k));
  }

  static float64_t
  _rd_exp(float64_t x)
  {
    union doub r0;
    //  degree-11 minimax coeffs c0..c11 (math.hoon:1542-1547)
    static const c3_d cs[12] = {
      0x3ff0000000000000ULL, 0x3ff0000000000000ULL, 0x3fe0000000000011ULL,
      0x3fc555555555555aULL, 0x3fa555555554f0cfULL, 0x3f8111111110f225ULL,
      0x3f56c16c187fbe02ULL, 0x3f2a01a01b14378fULL, 0x3efa01991ac8730aULL,
      0x3ec71ddf5749d126ULL, 0x3e928b4057f44145ULL, 0x3e5af631d0059becULL
    };
    union doub log2e, ln2hi, ln2lo, ka, kf, rr, p, c, zero;

    zero.c = 0;
    r0.d = x;
    if ( !f64_eq(x, x) )       { r0.c = _RD_QNAN; return r0.d; }  // NaN
    if ( r0.c == _RD_PINF )    { return x; }                      // +inf
    if ( r0.c == _RD_NINF )    { r0.c = 0; return r0.d; }         // -inf -> 0

    log2e.c = 0x3ff71547652b82feULL;
    ln2hi.c = 0x3fe62e42fee00000ULL;
    ln2lo.c = 0x3dea39ef35793c76ULL;

    c3_ds k = (c3_ds)f64_to_i64(f64_mul(x, log2e.d), softfloat_round_near_even, 0);
    if ( (k - 1025) >= 0 )    { r0.c = _RD_PINF; return r0.d; }   // overflow -> inf
    if ( !((k + 1075) >= 0) ) { r0.c = 0; return r0.d; }          // underflow -> 0

    ka.d = ui64_to_f64( (c3_d)(k < 0 ? -k : k) );
    kf.d = (k >= 0) ? ka.d : f64_sub(zero.d, ka.d);
    rr.d = f64_sub( f64_sub(x, f64_mul(kf.d, ln2hi.d)), f64_mul(kf.d, ln2lo.d) );

    p.c = 0;
    for ( int i = 12; i-- != 0; ) {        // Horner over flop(cs): c11..c0
      c.c = cs[i];
      p.d = f64_add(f64_mul(p.d, rr.d), c.d);
    }
    return _rd_scale2(p.d, k);
  }

/* @rd log/log-2/log-10 -- math.hoon ++rd ++log/++log-2/++log-10/++lr
**   lr: x = 2^e * m, m in [sqrt(1/2),sqrt(2)); returns [e-as-rd, log(1+f)] with
**   f=m-1, s=f/(2+f), log(1+f) = f - s*(f - 2z*P2(z)), z=s*s, P2 the atanh
**   series 1/3 + z/5 + ...  log/log-2/log-10 recombine the integer part.
*/
  //  +lr (math.hoon:2043): for finite positive x, -> *ef = e as @rd, *lm = log(m)
  static void _rd_lr(float64_t x, float64_t* ef, float64_t* lm) {
    static const c3_d cs[10] = {
      0x3fd5555555555555ULL, 0x3fc999999999999aULL, 0x3fc2492492492492ULL,
      0x3fbc71c71c71c71cULL, 0x3fb745d1745d1746ULL, 0x3fb3b13b13b13b14ULL,
      0x3fb1111111111111ULL, 0x3fae1e1e1e1e1e1eULL, 0x3faaf286bca1af28ULL,
      0x3fa8618618618618ULL };
    union doub r0, xx, m, f, s, z, p, c, rr, l1, efd, one, zero;
    zero.c = 0; one.c = 0x3ff0000000000000ULL;
    r0.d = x;
    int sub = ( ((r0.c >> 52) & 0x7ffULL) == 0 );
    xx = r0;
    if ( sub ) xx.d = f64_mul(x, _rd_bits(0x4350000000000000ULL).d);   // x * 2^54
    c3_ds ae = sub ? -54 : 0;
    c3_d  b  = xx.c;
    c3_ds e  = (c3_ds)((b >> 52) & 0x7ffULL) - 1023;
    m.c = (b & 0xfffffffffffffULL) | 0x3ff0000000000000ULL;            // m in [1,2)
    if ( f64_le(_rd_bits(0x3ff6a09e667f3bcdULL).d, m.d) ) {            // m >= sqrt(2)
      m.d = f64_mul(m.d, _rd_bits(0x3fe0000000000000ULL).d); e = e + 1;
    }
    e = e + ae;
    f.d = f64_sub(m.d, one.d);
    s.d = f64_div(f.d, f64_add(m.d, one.d));
    z.d = f64_mul(s.d, s.d);
    p.c = 0;
    for ( int i = 10; i-- != 0; ) { c.c = cs[i]; p.d = f64_add(f64_mul(p.d, z.d), c.d); }
    rr.d = f64_mul(f64_add(z.d, z.d), p.d);
    l1.d = f64_sub(f.d, f64_mul(s.d, f64_sub(f.d, rr.d)));
    efd.d = ui64_to_f64( (c3_d)(e < 0 ? -e : e) );
    if ( e < 0 ) efd.d = f64_sub(zero.d, efd.d);
    *ef = efd.d; *lm = l1.d;
  }
  //  shared guard: 0 on NaN/+inf/+-0/x<0 returned via *out
  static int _rd_log_guard(float64_t x, float64_t* out) {
    union doub r0; r0.d = x;
    if ( !f64_eq(x, x) )                       { r0.c = _RD_QNAN; *out = r0.d; return 1; }
    if ( r0.c == _RD_PINF )                    { *out = x;        return 1; }
    if ( (r0.c == 0)||(r0.c == 0x8000000000000000ULL) ) { r0.c = _RD_NINF; *out = r0.d; return 1; }
    if ( (r0.c >> 63) == 1 )                   { r0.c = _RD_QNAN; *out = r0.d; return 1; }
    return 0;
  }
  static float64_t _rd_log(float64_t x) {
    float64_t g, ef, lm; union doub hi, lo;
    if ( _rd_log_guard(x, &g) ) return g;
    _rd_lr(x, &ef, &lm);
    hi.d = f64_mul(ef, _rd_bits(0x3fe62e42fee00000ULL).d);            // e*ln2hi
    lo.d = f64_mul(ef, _rd_bits(0x3dea39ef35793c76ULL).d);            // e*ln2lo
    return f64_add(hi.d, f64_add(lm, lo.d));
  }
  static float64_t _rd_log2(float64_t x) {
    float64_t g, ef, lm;
    if ( _rd_log_guard(x, &g) ) return g;
    _rd_lr(x, &ef, &lm);                                              // e + lm/ln2
    return f64_add(ef, f64_mul(lm, _rd_bits(0x3ff71547652b82feULL).d));
  }
  static float64_t _rd_log10(float64_t x) {
    float64_t g, ef, lm;
    if ( _rd_log_guard(x, &g) ) return g;
    _rd_lr(x, &ef, &lm);                                              // e*log10(2) + lm/ln10
    return f64_add(f64_mul(ef, _rd_bits(0x3fd34413509f79ffULL).d),
                   f64_mul(lm, _rd_bits(0x3fdbcb7b1526e50eULL).d));
  }

  static inline float64_t _rd_neg(float64_t a) { union doub z; z.c=0; return f64_sub(z.d, a); }

/* @rd sin/cos -- math.hoon ++rd ++sin/++cos/++rd-trig
**   x = q*(pi/2) + (rhi+rlo) (2-part pi/2), fdlibm sin/cos kernels by q&3.
*/
  static const c3_d _RD_SC[8] = {     // sin kernel coeffs (math.hoon:1611)
    0xbfc5555555555555ULL, 0x3f81111111111111ULL, 0xbf2a01a01a01a01aULL,
    0x3ec71de3a556c734ULL, 0xbe5ae64567f544e4ULL, 0x3de6124613a86d09ULL,
    0xbd6ae7f3e733b81fULL, 0x3ce952c77030ad4aULL
  };
  static const c3_d _RD_CC[8] = {     // cos kernel coeffs (math.hoon:1618)
    0x3fa5555555555555ULL, 0xbf56c16c16c16c17ULL, 0x3efa01a01a01a01aULL,
    0xbe927e4fb7789f5cULL, 0x3e21eed8eff8d898ULL, 0xbda93974a8c07c9dULL,
    0x3d2ae7f3e733b81fULL, 0xbca6827863b97d97ULL
  };

  static float64_t _rd_ksin(float64_t xx, float64_t yy) {
    union doub z, r, v, aa, bb, dd, c, half;
    half.c = 0x3fe0000000000000ULL;
    z.d = f64_mul(xx, xx);
    r.c = 0;                            // Horner over flop(tail sc): sc[7..1]
    for ( int i = 8; i-- != 1; ) { c.c = _RD_SC[i]; r.d = f64_add(f64_mul(r.d, z.d), c.d); }
    v.d = f64_mul(z.d, xx);
    aa.d = f64_sub(f64_mul(half.d, yy), f64_mul(v.d, r.d));
    bb.d = f64_sub(f64_mul(z.d, aa.d), yy);
    dd.d = f64_sub(bb.d, f64_mul(v.d, _rd_bits(_RD_SC[0]).d));
    return f64_sub(xx, dd.d);
  }
  static float64_t _rd_kcos(float64_t xx, float64_t yy) {
    union doub z, rc, hz, w2, aa, bb, c, half, one;
    half.c = 0x3fe0000000000000ULL; one.c = 0x3ff0000000000000ULL;
    z.d = f64_mul(xx, xx);
    rc.c = 0;                          // Horner over flop(cc): cc[7..0]
    for ( int i = 8; i-- != 0; ) { c.c = _RD_CC[i]; rc.d = f64_add(f64_mul(rc.d, z.d), c.d); }
    hz.d = f64_mul(half.d, z.d);
    w2.d = f64_sub(one.d, hz.d);
    aa.d = f64_sub(f64_sub(one.d, w2.d), hz.d);
    bb.d = f64_sub(f64_mul(f64_mul(z.d, z.d), rc.d), f64_mul(xx, yy));
    return f64_add(w2.d, f64_add(aa.d, bb.d));
  }
  //  trig-fin: is_sin ? sin(x) : cos(x); ax=|x|, sb=sign bit (math.hoon:1643)
  static float64_t _rd_trigfin(int is_sin, float64_t ax, c3_d sb) {
    union doub qf, t, w, rhi, rlo, ks, kc, v;
    c3_ds q = (c3_ds)f64_to_i64(f64_mul(ax, _rd_bits(0x3fe45f306dc9c883ULL).d),
                                softfloat_round_near_even, 0);          // round(ax*2/pi)
    c3_d  aq = (c3_d)(q < 0 ? -q : q);
    qf.d = ui64_to_f64(aq);
    t.d = f64_sub(ax, f64_mul(qf.d, _rd_bits(0x3ff921fb54400000ULL).d)); // ax - qf*pi/2_hi
    w.d = f64_mul(qf.d, _rd_bits(0x3dd0b4611a626331ULL).d);             // qf*pi/2_lo
    rhi.d = f64_sub(t.d, w.d);
    rlo.d = f64_sub(f64_sub(t.d, rhi.d), w.d);
    int m = (int)(aq & 3);
    ks.d = _rd_ksin(rhi.d, rlo.d);
    kc.d = _rd_kcos(rhi.d, rlo.d);
    if ( is_sin ) {
      v.d = (m==0) ? ks.d : (m==1) ? kc.d : (m==2) ? _rd_neg(ks.d) : _rd_neg(kc.d);
      return (sb == 1) ? _rd_neg(v.d) : v.d;
    }
    return (m==0) ? kc.d : (m==1) ? _rd_neg(ks.d) : (m==2) ? _rd_neg(kc.d) : ks.d;
  }
  static float64_t _rd_sin(float64_t x) {
    union doub r0, ax;
    r0.d = x;
    if ( !f64_eq(x, x) )                   { r0.c = _RD_QNAN; return r0.d; }  // NaN
    if ( (r0.c == _RD_PINF)||(r0.c == _RD_NINF) ) { r0.c = _RD_QNAN; return r0.d; }  // +-inf -> NaN
    if ( (r0.c == 0)||(r0.c == 0x8000000000000000ULL) ) { return x; }         // +-0 -> +-0
    ax.c = r0.c & 0x7fffffffffffffffULL;
    return _rd_trigfin(1, ax.d, r0.c >> 63);
  }
  static float64_t _rd_cos(float64_t x) {
    union doub r0, ax;
    r0.d = x;
    if ( !f64_eq(x, x) )                   { r0.c = _RD_QNAN; return r0.d; }  // NaN
    if ( (r0.c == _RD_PINF)||(r0.c == _RD_NINF) ) { r0.c = _RD_QNAN; return r0.d; }  // +-inf -> NaN
    ax.c = r0.c & 0x7fffffffffffffffULL;
    return _rd_trigfin(0, ax.d, 0);
  }

/* @rd tan -- math.hoon ++rd ++tan/++rd-tan (fdlibm __kernel_tan)
**   q*pi/2 reduction; big |x|~pi/4 reduced; odd q -> -cot path.
*/
  static float64_t _rd_ktan(float64_t x, float64_t y, c3_ds iy) {
    static const c3_d rl[6] = {        // math.hoon:1700
      0x3fc111111110fe7aULL, 0x3f9664f48406d637ULL, 0x3f6d6d22c9560328ULL,
      0x3f4344d8f2f26501ULL, 0x3f147e88a03792a6ULL, 0xbef375cbdb605373ULL };
    static const c3_d vl[6] = {        // math.hoon:1705
      0x3faba1ba1bb341feULL, 0x3f8226e3e96e8493ULL, 0x3f57dbc8fee08315ULL,
      0x3f3026f71a8d1068ULL, 0x3f12b80f32f0a7e9ULL, 0x3efb2a7074bf7ad4ULL };
    union doub xb, ax, xa, ya, xr, yr, z, w, c, rr, vp, vv, s, r, w2;
    union doub one, mone, two, third;
    one.c=0x3ff0000000000000ULL; mone.c=0xbff0000000000000ULL;
    two.c=0x4000000000000000ULL;  third.c=0x3fd5555555555563ULL;
    xb.d = x;
    c3_d hxneg = xb.c >> 63;
    ax.c = xb.c & 0x7fffffffffffffffULL;
    int big = f64_le(_rd_bits(0x3fe5942800000000ULL).d, ax.d);        // |x| >= ~0.674
    xa.d = (hxneg == 1) ? _rd_neg(x) : x;
    ya.d = (hxneg == 1) ? _rd_neg(y) : y;
    if ( big ) {                       // pio4_hi - xa + (pio4_lo - ya)
      xr.d = f64_add(f64_sub(_rd_bits(0x3fe921fb54442d18ULL).d, xa.d),
                     f64_sub(_rd_bits(0x3c81a62633145c07ULL).d, ya.d));
      yr.c = 0;
    } else { xr.d = x; yr.d = y; }
    z.d = f64_mul(xr.d, xr.d);
    w.d = f64_mul(z.d, z.d);
    rr.c = 0;
    for ( int i = 6; i-- != 0; ) { c.c = rl[i]; rr.d = f64_add(f64_mul(rr.d, w.d), c.d); }
    vp.c = 0;
    for ( int i = 6; i-- != 0; ) { c.c = vl[i]; vp.d = f64_add(f64_mul(vp.d, w.d), c.d); }
    vv.d = f64_mul(z.d, vp.d);
    s.d = f64_mul(z.d, xr.d);
    r.d = f64_add(yr.d, f64_mul(z.d, f64_add(f64_mul(s.d, f64_add(rr.d, vv.d)), yr.d)));
    r.d = f64_add(r.d, f64_mul(third.d, s.d));
    w2.d = f64_add(xr.d, r.d);
    if ( big ) {
      union doub fac, v, t1;
      fac.d = (hxneg == 1) ? mone.d : one.d;
      v.d   = (iy == 1) ? one.d : mone.d;
      t1.d = f64_sub(f64_div(f64_mul(w2.d, w2.d), f64_add(w2.d, v.d)), r.d);
      t1.d = f64_mul(two.d, f64_sub(xr.d, t1.d));
      return f64_mul(fac.d, f64_sub(v.d, t1.d));
    }
    if ( iy == 1 ) return w2.d;
    { union doub zz, vv2, a, tt, ss;     // -cot path
      zz.c  = w2.c & 0xffffffff00000000ULL;
      vv2.d = f64_sub(r.d, f64_sub(zz.d, xr.d));
      a.d   = f64_div(mone.d, w2.d);
      tt.c  = a.c & 0xffffffff00000000ULL;
      ss.d  = f64_add(one.d, f64_mul(tt.d, zz.d));
      return f64_add(tt.d, f64_mul(a.d, f64_add(ss.d, f64_mul(tt.d, vv2.d))));
    }
  }
  static float64_t _rd_tan(float64_t x) {
    union doub r0, ax, qf, t, w, rhi, rlo, res;
    r0.d = x;
    if ( !f64_eq(x, x) )                   { r0.c = _RD_QNAN; return r0.d; }
    if ( (r0.c == _RD_PINF)||(r0.c == _RD_NINF) ) { r0.c = _RD_QNAN; return r0.d; }
    if ( (r0.c == 0)||(r0.c == 0x8000000000000000ULL) ) { return x; }
    c3_d neg = r0.c >> 63;
    ax.c = r0.c & 0x7fffffffffffffffULL;
    c3_ds q = (c3_ds)f64_to_i64(f64_mul(ax.d, _rd_bits(0x3fe45f306dc9c883ULL).d),
                                softfloat_round_near_even, 0);
    c3_d aq = (c3_d)(q < 0 ? -q : q);
    qf.d = ui64_to_f64(aq);
    t.d = f64_sub(ax.d, f64_mul(qf.d, _rd_bits(0x3ff921fb54400000ULL).d));
    w.d = f64_mul(qf.d, _rd_bits(0x3dd0b4611a626331ULL).d);
    rhi.d = f64_sub(t.d, w.d);
    rlo.d = f64_sub(f64_sub(t.d, rhi.d), w.d);
    c3_ds iy = ((aq & 1) == 0) ? 1 : -1;
    res.d = _rd_ktan(rhi.d, rlo.d, iy);
    return (neg == 1) ? _rd_neg(res.d) : res.d;
  }

/* @rd atan/atan2 -- math.hoon ++rd ++atan/++rd-atan/++atan2
**   fdlibm breakpoint reduction (7/16,11/16,19/16,39/16) + degree-10 minimax.
*/
  static float64_t _rd_atan(float64_t x) {
    static const c3_d at[11] = {       // fdlibm aT[] (math.hoon:1869)
      0x3fd555555555550dULL, 0xbfc999999998ebc4ULL, 0x3fc24924920083ffULL,
      0xbfbc71c6fe231671ULL, 0x3fb745cdc54c206eULL, 0xbfb3b0f2af749a6dULL,
      0x3fb10d66a0d03d51ULL, 0xbfadde2d52defd9aULL, 0x3fa97b4b24760debULL,
      0xbfa2b4442c6a6c2fULL, 0x3f90ad3ae322da11ULL };
    union doub r0, ax, xr, hi, lo, z, sp, s, c, res, one, two, ohf;
    r0.d = x;
    if ( !f64_eq(x, x) )       { r0.c = _RD_QNAN; return r0.d; }       // NaN
    if ( r0.c == _RD_PINF )    { r0.c = 0x3ff921fb54442d18ULL; return r0.d; }  // +inf -> pi/2
    if ( r0.c == _RD_NINF )    { r0.c = 0xbff921fb54442d18ULL; return r0.d; }  // -inf -> -pi/2
    if ( (r0.c == 0)||(r0.c == 0x8000000000000000ULL) ) { return x; }  // +-0
    c3_d neg = r0.c >> 63;
    ax.c = r0.c & 0x7fffffffffffffffULL;
    one.c=0x3ff0000000000000ULL; two.c=0x4000000000000000ULL; ohf.c=0x3ff8000000000000ULL;
    int dir = 0;
    if ( f64_lt(ax.d, _rd_bits(0x3fdc000000000000ULL).d) ) {            // |x| < 7/16
      xr.d = ax.d; hi.c = 0; lo.c = 0; dir = 1;
    } else if ( f64_lt(ax.d, _rd_bits(0x3fe6000000000000ULL).d) ) {     // < 11/16
      xr.d = f64_div(f64_sub(f64_add(ax.d, ax.d), one.d), f64_add(two.d, ax.d));
      hi.c = 0x3fddac670561bb4fULL; lo.c = 0x3c7a2b7f222f65e2ULL;       // atan(0.5)
    } else if ( f64_lt(ax.d, _rd_bits(0x3ff3000000000000ULL).d) ) {     // < 19/16
      xr.d = f64_div(f64_sub(ax.d, one.d), f64_add(ax.d, one.d));
      hi.c = 0x3fe921fb54442d18ULL; lo.c = 0x3c81a62633145c07ULL;       // atan(1)=pi/4
    } else if ( f64_lt(ax.d, _rd_bits(0x4003800000000000ULL).d) ) {     // < 39/16
      xr.d = f64_div(f64_sub(ax.d, ohf.d), f64_add(one.d, f64_mul(ohf.d, ax.d)));
      hi.c = 0x3fef730bd281f69bULL; lo.c = 0x3c7007887af0cbbdULL;       // atan(1.5)
    } else {                                                            // -1/x
      xr.d = f64_div(_rd_bits(0xbff0000000000000ULL).d, ax.d);
      hi.c = 0x3ff921fb54442d18ULL; lo.c = 0x3c91a62633145c07ULL;       // pi/2
    }
    z.d = f64_mul(xr.d, xr.d);
    sp.c = 0;
    for ( int i = 11; i-- != 0; ) { c.c = at[i]; sp.d = f64_add(f64_mul(sp.d, z.d), c.d); }
    s.d = f64_mul(z.d, sp.d);
    if ( dir ) res.d = f64_sub(xr.d, f64_mul(xr.d, s.d));
    else       res.d = f64_sub(hi.d, f64_sub(f64_sub(f64_mul(xr.d, s.d), lo.d), xr.d));
    return (neg == 1) ? _rd_neg(res.d) : res.d;
  }
  static float64_t _rd_atan2(float64_t y, float64_t x) {
    union doub xb, pi, pi2, mpi2, zero;
    zero.c = 0; pi.c = 0x400921fb54442d18ULL;
    pi2.c = 0x3ff921fb54442d18ULL; mpi2.c = 0xbff921fb54442d18ULL;
    xb.d = x;
    if ( f64_lt(zero.d, x) ) return _rd_atan(f64_div(y, x));            // x>0
    if ( f64_lt(x, zero.d) && f64_le(zero.d, y) )                       // x<0,y>=0
      return f64_add(_rd_atan(f64_div(y, x)), pi.d);
    if ( f64_lt(x, zero.d) && f64_lt(y, zero.d) )                       // x<0,y<0
      return f64_sub(_rd_atan(f64_div(y, x)), pi.d);
    if ( (xb.c == 0) && f64_lt(zero.d, y) ) return pi2.d;               // x==+0,y>0
    if ( (xb.c == 0) && f64_lt(y, zero.d) ) return mpi2.d;              // x==+0,y<0
    return zero.d;
  }

/* @rd asin/acos -- math.hoon ++rd ++asin/++acos/++rd-ainv
**   fdlibm rational P/Q kernel; sqt is correctly-rounded (= f64_sqrt).
*/
  static float64_t _rd_ainv_rr(float64_t t) {   // R(t) = P(t)/Q(t)  (math.hoon:1775)
    static const c3_d ps[6] = {
      0x3fc5555555555555ULL, 0xbfd4d61203eb6f7dULL, 0x3fc9c1550e884455ULL,
      0xbfa48228b5688f3bULL, 0x3f49efe07501b288ULL, 0x3f023de10dfdf709ULL };
    static const c3_d qs[4] = {
      0xc0033a271c8a2d4bULL, 0x40002ae59c598ac8ULL, 0xbfe6066c1b8d0159ULL,
      0x3fb3b8c5b12e9282ULL };
    union doub pp, qq, c, one;
    one.c = 0x3ff0000000000000ULL;
    pp.c = 0;
    for ( int i = 6; i-- != 0; ) { c.c = ps[i]; pp.d = f64_add(f64_mul(pp.d, t), c.d); }
    qq.c = 0;
    for ( int i = 4; i-- != 0; ) { c.c = qs[i]; qq.d = f64_add(f64_mul(qq.d, t), c.d); }
    return f64_div(f64_mul(t, pp.d), f64_add(one.d, f64_mul(t, qq.d)));
  }
  static float64_t _rd_asin(float64_t x) {
    union doub r0, ax, t, w, r, s, res, half, one, two, pio2h, pio2l, pio4;
    half.c=0x3fe0000000000000ULL; one.c=0x3ff0000000000000ULL; two.c=0x4000000000000000ULL;
    pio2h.c=0x3ff921fb54442d18ULL; pio2l.c=0x3c91a62633145c07ULL; pio4.c=0x3fe921fb54442d18ULL;
    r0.d = x;
    if ( !f64_eq(x, x) )            { r0.c = _RD_QNAN; return r0.d; }   // NaN
    c3_d sgn = r0.c >> 63;
    ax.c = r0.c & 0x7fffffffffffffffULL;
    if ( f64_lt(one.d, ax.d) )      { r0.c = _RD_QNAN; return r0.d; }   // |x|>1 -> NaN
    if ( ax.c == one.c )                                               // |x|==1
      return f64_add(f64_mul(x, pio2h.d), f64_mul(x, pio2l.d));
    if ( f64_lt(ax.d, half.d) ) {                                     // |x|<0.5
      if ( f64_lt(ax.d, _rd_bits(0x3e50000000000000ULL).d) ) return x; // tiny
      t.d = f64_mul(x, x);
      return f64_add(x, f64_mul(x, _rd_ainv_rr(t.d)));
    }
    w.d = f64_sub(one.d, ax.d);
    t.d = f64_mul(w.d, half.d);
    r.d = _rd_ainv_rr(t.d);
    s.d = f64_sqrt(t.d);
    if ( f64_le(_rd_bits(0x3fef333300000000ULL).d, ax.d) ) {           // near 1
      res.d = f64_sub(pio2h.d, f64_sub(f64_mul(two.d, f64_add(s.d, f64_mul(s.d, r.d))), pio2l.d));
      return (sgn == 1) ? _rd_neg(res.d) : res.d;
    }
    { union doub df, cc, p2, q2;                                       // head/tail
      df.c = s.c & 0xffffffff00000000ULL;
      cc.d = f64_div(f64_sub(t.d, f64_mul(df.d, df.d)), f64_add(s.d, df.d));
      p2.d = f64_sub(f64_mul(two.d, f64_mul(s.d, r.d)), f64_sub(pio2l.d, f64_mul(two.d, cc.d)));
      q2.d = f64_sub(pio4.d, f64_mul(two.d, df.d));
      res.d = f64_sub(pio4.d, f64_sub(p2.d, q2.d));
      return (sgn == 1) ? _rd_neg(res.d) : res.d;
    }
  }
  static float64_t _rd_acos(float64_t x) {
    union doub r0, ax, z, s, r, w, half, one, two, pi, pio2h, pio2l;
    half.c=0x3fe0000000000000ULL; one.c=0x3ff0000000000000ULL; two.c=0x4000000000000000ULL;
    pi.c=0x400921fb54442d18ULL; pio2h.c=0x3ff921fb54442d18ULL; pio2l.c=0x3c91a62633145c07ULL;
    r0.d = x;
    if ( !f64_eq(x, x) )            { r0.c = _RD_QNAN; return r0.d; }   // NaN
    c3_d neg = r0.c >> 63;
    ax.c = r0.c & 0x7fffffffffffffffULL;
    if ( f64_lt(one.d, ax.d) )      { r0.c = _RD_QNAN; return r0.d; }   // |x|>1 -> NaN
    if ( ax.c == one.c ) {                                             // |x|==1
      if ( neg == 0 ) { union doub z0; z0.c = 0; return z0.d; }         // 1 -> 0
      return f64_add(pi.d, f64_mul(two.d, pio2l.d));                    // -1 -> pi
    }
    if ( f64_lt(ax.d, half.d) ) {                                     // |x|<0.5
      if ( f64_lt(ax.d, _rd_bits(0x3c60000000000000ULL).d) ) return pio2h.d;  // tiny -> pi/2
      z.d = f64_mul(x, x);
      r.d = _rd_ainv_rr(z.d);
      return f64_sub(pio2h.d, f64_sub(x, f64_sub(pio2l.d, f64_mul(x, r.d))));
    }
    if ( neg == 1 ) {                                                 // x <= -0.5
      z.d = f64_mul(f64_add(one.d, x), half.d);
      s.d = f64_sqrt(z.d);
      r.d = _rd_ainv_rr(z.d);
      w.d = f64_sub(f64_mul(r.d, s.d), pio2l.d);
      return f64_sub(pi.d, f64_mul(two.d, f64_add(s.d, w.d)));
    }
    { union doub df, cc;                                              // x >= 0.5
      z.d = f64_mul(f64_sub(one.d, x), half.d);
      s.d = f64_sqrt(z.d);
      df.c = s.c & 0xffffffff00000000ULL;
      cc.d = f64_div(f64_sub(z.d, f64_mul(df.d, df.d)), f64_add(s.d, df.d));
      r.d = _rd_ainv_rr(z.d);
      w.d = f64_add(f64_mul(r.d, s.d), cc.d);
      return f64_mul(two.d, f64_add(df.d, w.d));
    }
  }

/* @rd sqt/cbt -- math.hoon ++rd ++sqt/++cbt
**   sqt: correctly-rounded f64 sqrt (the Markstein-corrected Hoon == f64_sqrt).
**   cbt: sign(x) * exp(log|x| / 3).
*/
  static float64_t _rd_sqt(float64_t x) {
    union doub r0; r0.d = x;
    if ( !f64_eq(x, x) )                       { r0.c = _RD_QNAN; return r0.d; }  // NaN
    if ( r0.c == _RD_PINF )                    { return x; }                      // +inf
    if ( (r0.c == 0)||(r0.c == 0x8000000000000000ULL) ) { return x; }             // +-0
    if ( (r0.c >> 63) == 1 )                   { r0.c = _RD_QNAN; return r0.d; }  // x<0 -> NaN
    return f64_sqrt(x);
  }
  static float64_t _rd_cbt(float64_t x) {
    union doub r0, ax, r;
    r0.d = x;
    if ( !f64_eq(x, x) )                       { return x; }                      // NaN -> NaN
    if ( (r0.c == 0)||(r0.c == 0x8000000000000000ULL) ) { return x; }             // +-0
    ax.c = r0.c & 0x7fffffffffffffffULL;
    r.d = _rd_exp(f64_mul(_rd_log(ax.d), _rd_bits(0x3fd5555555555555ULL).d));     // exp(log|x|/3)
    return ((r0.c >> 63) == 1) ? _rd_neg(r.d) : r.d;
  }

/* @rd pow/pow-n -- math.hoon ++rd ++pow/++pow-n
**   pow-n: x^n by repeated mul (n a positive integer as @rd).
**   pow: positive-integer fast path -> pow-n, else exp(n*log x).
*/
  static float64_t _rd_pow_n(float64_t x, float64_t n) {
    union doub nn, p, one, two;
    one.c = 0x3ff0000000000000ULL; two.c = 0x4000000000000000ULL;
    nn.d = n;
    if ( nn.c == 0 ) return one.d;             // n == +0 -> 1
    p.d = x;
    while ( !f64_lt(n, two.d) ) {              // while n >= 2: p *= x; n -= 1
      p.d = f64_mul(p.d, x);
      n = f64_sub(n, one.d);
    }
    return p.d;
  }
  static float64_t _rd_pow(float64_t x, float64_t n) {
    union doub nn, ni, zero;
    zero.c = 0; nn.d = n;
    ni.d = i64_to_f64(f64_to_i64(n, softfloat_round_near_even, 0));   // san (need (toi n))
    if ( (nn.c == ni.c) && f64_lt(zero.d, n) )                        // positive integer
      return _rd_pow_n(x, ni.d);
    return _rd_exp(f64_mul(n, _rd_log(x)));                           // exp(n*log x)
  }

#ifndef MATH_JET_HARNESS

/* u3 ABI wrappers.  Each transcendental forces round-near-even; the @rd door's
** rounding axis is ignored (the Hoon does the same).  Sample is the single @rd
** at u3x_sam.
*/
  static u3_noun
  _rd_jet(u3_noun cor, float64_t (*fun)(float64_t))
  {
    u3_noun x = u3r_at(u3x_sam, cor);
    if ( u3_none == x || c3n == u3ud(x) ) {
      return u3m_bail(c3__exit);
    }
    {
      union doub c, e;
      softfloat_roundingMode = softfloat_round_near_even;
      c.c = u3r_chub(0, x);
      e.d = fun(c.d);
      return u3i_chubs(1, &e.c);
    }
  }

  u3_noun u3qi_rd_exp(u3_atom a)
  {
    union doub c, e;
    softfloat_roundingMode = softfloat_round_near_even;
    c.c = u3r_chub(0, a);
    e.d = _rd_exp(c.d);
    return u3i_chubs(1, &e.c);
  }
  u3_noun u3wi_rd_exp(u3_noun cor) { return _rd_jet(cor, _rd_exp); }

  u3_noun u3qi_rd_log(u3_atom a)
  {
    union doub c, e;
    softfloat_roundingMode = softfloat_round_near_even;
    c.c = u3r_chub(0, a);
    e.d = _rd_log(c.d);
    return u3i_chubs(1, &e.c);
  }
  u3_noun u3wi_rd_log(u3_noun cor) { return _rd_jet(cor, _rd_log); }

  u3_noun u3qi_rd_sin(u3_atom a)
  {
    union doub c, e;
    softfloat_roundingMode = softfloat_round_near_even;
    c.c = u3r_chub(0, a);
    e.d = _rd_sin(c.d);
    return u3i_chubs(1, &e.c);
  }
  u3_noun u3wi_rd_sin(u3_noun cor) { return _rd_jet(cor, _rd_sin); }

  u3_noun u3qi_rd_cos(u3_atom a)
  {
    union doub c, e;
    softfloat_roundingMode = softfloat_round_near_even;
    c.c = u3r_chub(0, a);
    e.d = _rd_cos(c.d);
    return u3i_chubs(1, &e.c);
  }
  u3_noun u3wi_rd_cos(u3_noun cor) { return _rd_jet(cor, _rd_cos); }

  u3_noun u3qi_rd_tan(u3_atom a)
  {
    union doub c, e;
    softfloat_roundingMode = softfloat_round_near_even;
    c.c = u3r_chub(0, a);
    e.d = _rd_tan(c.d);
    return u3i_chubs(1, &e.c);
  }
  u3_noun u3wi_rd_tan(u3_noun cor) { return _rd_jet(cor, _rd_tan); }

  u3_noun u3qi_rd_atan(u3_atom a)
  {
    union doub c, e;
    softfloat_roundingMode = softfloat_round_near_even;
    c.c = u3r_chub(0, a);
    e.d = _rd_atan(c.d);
    return u3i_chubs(1, &e.c);
  }
  u3_noun u3wi_rd_atan(u3_noun cor) { return _rd_jet(cor, _rd_atan); }

  u3_noun u3qi_rd_atan2(u3_atom y, u3_atom x)
  {
    union doub yy, xx, e;
    softfloat_roundingMode = softfloat_round_near_even;
    yy.c = u3r_chub(0, y);
    xx.c = u3r_chub(0, x);
    e.d = _rd_atan2(yy.d, xx.d);
    return u3i_chubs(1, &e.c);
  }
  u3_noun u3wi_rd_atan2(u3_noun cor)
  {
    u3_noun y, x;
    if ( c3n == u3r_mean(cor, {u3x_sam_2, &y}, {u3x_sam_3, &x}) ||
         c3n == u3ud(y) || c3n == u3ud(x) ) {
      return u3m_bail(c3__exit);
    }
    return u3qi_rd_atan2(y, x);
  }

  u3_noun u3qi_rd_asin(u3_atom a)
  {
    union doub c, e;
    softfloat_roundingMode = softfloat_round_near_even;
    c.c = u3r_chub(0, a);
    e.d = _rd_asin(c.d);
    return u3i_chubs(1, &e.c);
  }
  u3_noun u3wi_rd_asin(u3_noun cor) { return _rd_jet(cor, _rd_asin); }

  u3_noun u3qi_rd_acos(u3_atom a)
  {
    union doub c, e;
    softfloat_roundingMode = softfloat_round_near_even;
    c.c = u3r_chub(0, a);
    e.d = _rd_acos(c.d);
    return u3i_chubs(1, &e.c);
  }
  u3_noun u3wi_rd_acos(u3_noun cor) { return _rd_jet(cor, _rd_acos); }

  u3_noun u3qi_rd_log2(u3_atom a)
  { union doub c,e; softfloat_roundingMode=softfloat_round_near_even;
    c.c=u3r_chub(0,a); e.d=_rd_log2(c.d); return u3i_chubs(1,&e.c); }
  u3_noun u3wi_rd_log2(u3_noun cor) { return _rd_jet(cor, _rd_log2); }

  u3_noun u3qi_rd_log10(u3_atom a)
  { union doub c,e; softfloat_roundingMode=softfloat_round_near_even;
    c.c=u3r_chub(0,a); e.d=_rd_log10(c.d); return u3i_chubs(1,&e.c); }
  u3_noun u3wi_rd_log10(u3_noun cor) { return _rd_jet(cor, _rd_log10); }

  u3_noun u3qi_rd_sqt(u3_atom a)
  { union doub c,e; softfloat_roundingMode=softfloat_round_near_even;
    c.c=u3r_chub(0,a); e.d=_rd_sqt(c.d); return u3i_chubs(1,&e.c); }
  u3_noun u3wi_rd_sqt(u3_noun cor) { return _rd_jet(cor, _rd_sqt); }

  u3_noun u3qi_rd_cbt(u3_atom a)
  { union doub c,e; softfloat_roundingMode=softfloat_round_near_even;
    c.c=u3r_chub(0,a); e.d=_rd_cbt(c.d); return u3i_chubs(1,&e.c); }
  u3_noun u3wi_rd_cbt(u3_noun cor) { return _rd_jet(cor, _rd_cbt); }

  //  pow / pow-n: sample is [x=@rd n=@rd]
  static u3_noun _rd_jet2(u3_noun cor, float64_t (*fun)(float64_t, float64_t))
  {
    u3_noun x, n;
    if ( c3n == u3r_mean(cor, {u3x_sam_2, &x}, {u3x_sam_3, &n}) ||
         c3n == u3ud(x) || c3n == u3ud(n) ) {
      return u3m_bail(c3__exit);
    }
    { union doub xx, nn, e;
      softfloat_roundingMode = softfloat_round_near_even;
      xx.c = u3r_chub(0, x); nn.c = u3r_chub(0, n);
      e.d = fun(xx.d, nn.d);
      return u3i_chubs(1, &e.c);
    }
  }
  u3_noun u3qi_rd_pow(u3_atom x, u3_atom n)
  { union doub xx,nn,e; softfloat_roundingMode=softfloat_round_near_even;
    xx.c=u3r_chub(0,x); nn.c=u3r_chub(0,n); e.d=_rd_pow(xx.d,nn.d); return u3i_chubs(1,&e.c); }
  u3_noun u3wi_rd_pow(u3_noun cor) { return _rd_jet2(cor, _rd_pow); }

  u3_noun u3qi_rd_pow_n(u3_atom x, u3_atom n)
  { union doub xx,nn,e; softfloat_roundingMode=softfloat_round_near_even;
    xx.c=u3r_chub(0,x); nn.c=u3r_chub(0,n); e.d=_rd_pow_n(xx.d,nn.d); return u3i_chubs(1,&e.c); }
  u3_noun u3wi_rd_pow_n(u3_noun cor) { return _rd_jet2(cor, _rd_pow_n); }

#endif
