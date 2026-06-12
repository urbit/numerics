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

  //  The math doors carry a rounding mode r=?(%n %u %d %z) whose bunt is %z.
  //  The transcendental KERNELS force %n (correctly-rounded, no axis), but the
  //  composite arms (pow/atan2/tan/pow-n) round their BARE door ops per r.  The
  //  wrapper sets _math_rnd from the door's r (axis 60 of the gate); cores
  //  bracket bare ops with it and restore near-even around kernel calls.
  //  Default = minMag (%z), matching the door bunt for the harness/cold path.
  static int _math_rnd = softfloat_round_minMag;
  static inline int _rnd_of(c3_d r) {            // @tas 'n'/'u'/'d'/'z' -> SoftFloat
    switch ( r ) {
      case 'n': return softfloat_round_near_even;
      case 'u': return softfloat_round_max;
      case 'd': return softfloat_round_min;
      case 'z': return softfloat_round_minMag;
      default:  return softfloat_round_minMag;
    }
  }

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
  //  bare door ops (div/add/sub/mul) round per _math_rnd; atan kernel is %n.
  static float64_t _rd_atan2(float64_t y, float64_t x) {
    union doub xb, pi, two, mone, zero, q, a, r;
    zero.c = 0; pi.c = 0x400921fb54442d18ULL;
    two.c = 0x4000000000000000ULL; mone.c = 0xbff0000000000000ULL;
    xb.d = x;
    if ( f64_lt(zero.d, x) ) {                                          // x>0: atan(div y x)
      softfloat_roundingMode = _math_rnd; q.d = f64_div(y, x);
      softfloat_roundingMode = softfloat_round_near_even; return _rd_atan(q.d);
    }
    if ( f64_lt(x, zero.d) && f64_le(zero.d, y) ) {                     // x<0,y>=0: add(atan,pi)
      softfloat_roundingMode = _math_rnd; q.d = f64_div(y, x);
      softfloat_roundingMode = softfloat_round_near_even; a.d = _rd_atan(q.d);
      softfloat_roundingMode = _math_rnd; r.d = f64_add(a.d, pi.d);
      softfloat_roundingMode = softfloat_round_near_even; return r.d;
    }
    if ( f64_lt(x, zero.d) && f64_lt(y, zero.d) ) {                     // x<0,y<0: sub(atan,pi)
      softfloat_roundingMode = _math_rnd; q.d = f64_div(y, x);
      softfloat_roundingMode = softfloat_round_near_even; a.d = _rd_atan(q.d);
      softfloat_roundingMode = _math_rnd; r.d = f64_sub(a.d, pi.d);
      softfloat_roundingMode = softfloat_round_near_even; return r.d;
    }
    if ( (xb.c == 0) && f64_lt(zero.d, y) ) {                           // x==+0,y>0: div(pi,2)
      softfloat_roundingMode = _math_rnd; r.d = f64_div(pi.d, two.d);
      softfloat_roundingMode = softfloat_round_near_even; return r.d;
    }
    if ( (xb.c == 0) && f64_lt(y, zero.d) ) {                           // x==+0,y<0: mul(-1,div(pi,2))
      softfloat_roundingMode = _math_rnd;
      r.d = f64_mul(mone.d, f64_div(pi.d, two.d));
      softfloat_roundingMode = softfloat_round_near_even; return r.d;
    }
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
    softfloat_roundingMode = _math_rnd;        // bare mul/sub round per door r
    p.d = x;
    while ( !f64_lt(n, two.d) ) {              // while n >= 2: p *= x; n -= 1
      p.d = f64_mul(p.d, x);
      n = f64_sub(n, one.d);
    }
    softfloat_roundingMode = softfloat_round_near_even;
    return p.d;
  }
  static float64_t _rd_pow(float64_t x, float64_t n) {
    union doub nn, ni, zero, lg, prod;
    zero.c = 0; nn.d = n;
    //  integer detection is rounding-mode-independent (exact for true integers)
    ni.d = i64_to_f64(f64_to_i64(n, softfloat_round_near_even, 0));   // san (need (toi n))
    if ( (nn.c == ni.c) && f64_lt(zero.d, n) )                        // positive integer
      return _rd_pow_n(x, ni.d);
    lg.d = _rd_log(x);                                                // %n kernel
    softfloat_roundingMode = _math_rnd;                              // bare mul per door r
    prod.d = f64_mul(n, lg.d);
    softfloat_roundingMode = softfloat_round_near_even;
    return _rd_exp(prod.d);                                           // %n kernel: exp(n*log x)
  }

/* ===================================================================
** @rs (single-precision) cores -- math.hoon ++rs.  Single-precision twins
** of the @rd cores: identical reductions/Horner order, single-precision
** constants and minimax coeffs (read from the Hoon arms), SoftFloat f32 ops.
** Marshalling is still chub-based (read low 32 bits, write a 32-bit atom), so
** these are word-size-agnostic exactly like the @rd cores.  The bit pattern is
** a plain uint32_t (NOT c3_w, which is 64-bit on the vere64 build).
** =================================================================== */

  union sing {
    float32_t s;
    uint32_t  c;
  };

  static const uint32_t _RS_QNAN = 0x7fc00000U;
  static const uint32_t _RS_PINF = 0x7f800000U;
  static const uint32_t _RS_NINF = 0xff800000U;

  static inline union sing _rs_bits(uint32_t b) { union sing u; u.c = b; return u; }
  static inline float32_t  _rs_neg(float32_t a) {            // (sub .0 a)
    union sing z; z.c = 0; return f32_sub(z.s, a);
  }

/* @rs exp -- math.hoon ++rs ++exp
**   x = k*ln2 + r (Cody-Waite), exp(x) = 2^k * P(r), P a degree-6 minimax.
*/
  //  pow2(j) = (j+127)<<23 as f32 bits = 2^j  (normal range j in [-126,127])
  static inline float32_t _rs_pow2(c3_ds j) {
    union sing u; u.c = ((uint32_t)(j + 127)) << 23; return u.s;
  }
  //  scale2: ldexp with overflow/subnormal tails (math.hoon ++rs ++exp)
  static inline float32_t _rs_scale2(float32_t p, c3_ds k) {
    if ( (k - 128) >= 0 ) {                            // k>=128
      return f32_mul(f32_mul(p, _rs_pow2(127)), _rs_pow2(k - 127));
    }
    if ( !((k + 126) >= 0) ) {                         // k<-126
      return f32_mul(f32_mul(p, _rs_pow2(k + 24)), _rs_pow2(-24));
    }
    return f32_mul(p, _rs_pow2(k));
  }
  static float32_t _rs_exp(float32_t x) {
    union sing r0;
    //  degree-6 minimax coeffs c0..c6 (math.hoon ++rs ++exp)
    static const uint32_t cs[7] = {
      0x3f800000U, 0x3f800000U, 0x3f000000U, 0x3e2aaa02U,
      0x3d2aaa56U, 0x3c0937d3U, 0x3ab6ba99U
    };
    union sing log2e, ln2hi, ln2lo, ka, kf, rr, p, c, zero;
    zero.c = 0;
    r0.s = x;
    if ( !f32_eq(x, x) )    { r0.c = _RS_QNAN; return r0.s; }  // NaN
    if ( r0.c == _RS_PINF ) { return x; }                      // +inf
    if ( r0.c == _RS_NINF ) { r0.c = 0; return r0.s; }         // -inf -> 0

    log2e.c = 0x3fb8aa3bU;
    ln2hi.c = 0x3f317200U;
    ln2lo.c = 0x35bfbe8eU;

    c3_ds k = (c3_ds)f32_to_i32(f32_mul(x, log2e.s), softfloat_round_near_even, 0);
    if ( (k - 129) >= 0 )    { r0.c = _RS_PINF; return r0.s; }  // overflow -> inf
    if ( !((k + 150) >= 0) ) { r0.c = 0; return r0.s; }         // underflow -> 0

    ka.s = ui32_to_f32( (uint32_t)(k < 0 ? -k : k) );
    kf.s = (k >= 0) ? ka.s : f32_sub(zero.s, ka.s);
    rr.s = f32_sub( f32_sub(x, f32_mul(kf.s, ln2hi.s)), f32_mul(kf.s, ln2lo.s) );

    p.c = 0;
    for ( int i = 7; i-- != 0; ) {        // Horner over flop(cs): c6..c0
      c.c = cs[i];
      p.s = f32_add(f32_mul(p.s, rr.s), c.s);
    }
    return _rs_scale2(p.s, k);
  }

/* @rs sin/cos/tan -- math.hoon ++rs ++sin/++cos/++rs-trig/++tan
**   x = q*(pi/2) + (rhi+rlo) (3-part pi/2, f32 needs the extra word), fdlibm
**   sin/cos kernels by q&3.  tan = sin/cos (no dedicated kernel, unlike @rd).
*/
  static const uint32_t _RS_SC[5] = {     // sin kernel coeffs
    0xbe2aaaabU, 0x3c088889U, 0xb9500d01U, 0x3638ef1dU, 0xb2d7322bU
  };
  static const uint32_t _RS_CC[5] = {     // cos kernel coeffs
    0x3d2aaaabU, 0xbab60b61U, 0x37d00d01U, 0xb493f27eU, 0x310f76c7U
  };
  static float32_t _rs_ksin(float32_t xx, float32_t yy) {
    union sing z, r, v, aa, bb, dd, c, half;
    half.c = 0x3f000000U;
    z.s = f32_mul(xx, xx);
    r.c = 0;                            // Horner over flop(tail sc): sc[4..1]
    for ( int i = 5; i-- != 1; ) { c.c = _RS_SC[i]; r.s = f32_add(f32_mul(r.s, z.s), c.s); }
    v.s = f32_mul(z.s, xx);
    aa.s = f32_sub(f32_mul(half.s, yy), f32_mul(v.s, r.s));
    bb.s = f32_sub(f32_mul(z.s, aa.s), yy);
    dd.s = f32_sub(bb.s, f32_mul(v.s, _rs_bits(_RS_SC[0]).s));
    return f32_sub(xx, dd.s);
  }
  static float32_t _rs_kcos(float32_t xx, float32_t yy) {
    union sing z, rc, hz, w2, aa, bb, c, half, one;
    half.c = 0x3f000000U; one.c = 0x3f800000U;
    z.s = f32_mul(xx, xx);
    rc.c = 0;                          // Horner over flop(cc): cc[4..0]
    for ( int i = 5; i-- != 0; ) { c.c = _RS_CC[i]; rc.s = f32_add(f32_mul(rc.s, z.s), c.s); }
    hz.s = f32_mul(half.s, z.s);
    w2.s = f32_sub(one.s, hz.s);
    aa.s = f32_sub(f32_sub(one.s, w2.s), hz.s);
    bb.s = f32_sub(f32_mul(f32_mul(z.s, z.s), rc.s), f32_mul(xx, yy));
    return f32_add(w2.s, f32_add(aa.s, bb.s));
  }
  //  trig-fin: is_sin ? sin(x) : cos(x); ax=|x|, sb=sign bit
  static float32_t _rs_trigfin(int is_sin, float32_t ax, uint32_t sb) {
    union sing qf, r1, r2, w, rhi, rlo, ks, kc, v;
    c3_ds q = (c3_ds)f32_to_i32(f32_mul(ax, _rs_bits(0x3f22f983U).s),
                                softfloat_round_near_even, 0);          // round(ax*2/pi)
    c3_d  aq = (c3_d)(q < 0 ? -q : q);
    qf.s = ui32_to_f32((uint32_t)aq);
    r1.s = f32_sub(ax, f32_mul(qf.s, _rs_bits(0x3fc90000U).s));         // ax - qf*pio2_1
    r2.s = f32_sub(r1.s, f32_mul(qf.s, _rs_bits(0x39fda000U).s));       // r1 - qf*pio2_2
    w.s = f32_mul(qf.s, _rs_bits(0x33a22169U).s);                      // qf*pio2_3
    rhi.s = f32_sub(r2.s, w.s);
    rlo.s = f32_sub(f32_sub(r2.s, rhi.s), w.s);
    int m = (int)(aq & 3);
    ks.s = _rs_ksin(rhi.s, rlo.s);
    kc.s = _rs_kcos(rhi.s, rlo.s);
    if ( is_sin ) {
      v.s = (m==0) ? ks.s : (m==1) ? kc.s : (m==2) ? _rs_neg(ks.s) : _rs_neg(kc.s);
      return (sb == 1) ? _rs_neg(v.s) : v.s;
    }
    return (m==0) ? kc.s : (m==1) ? _rs_neg(ks.s) : (m==2) ? _rs_neg(kc.s) : ks.s;
  }
  static float32_t _rs_sin(float32_t x) {
    union sing r0, ax;
    r0.s = x;
    if ( !f32_eq(x, x) )                          { r0.c = _RS_QNAN; return r0.s; }  // NaN
    if ( (r0.c == _RS_PINF)||(r0.c == _RS_NINF) ) { r0.c = _RS_QNAN; return r0.s; }  // +-inf -> NaN
    if ( (r0.c == 0)||(r0.c == 0x80000000U) )     { return x; }                      // +-0 -> +-0
    ax.c = r0.c & 0x7fffffffU;
    return _rs_trigfin(1, ax.s, r0.c >> 31);
  }
  static float32_t _rs_cos(float32_t x) {
    union sing r0, ax;
    r0.s = x;
    if ( !f32_eq(x, x) )                          { r0.c = _RS_QNAN; return r0.s; }  // NaN
    if ( (r0.c == _RS_PINF)||(r0.c == _RS_NINF) ) { r0.c = _RS_QNAN; return r0.s; }  // +-inf -> NaN
    ax.c = r0.c & 0x7fffffffU;
    return _rs_trigfin(0, ax.s, 0);
  }
  //  tan = (div (sin x) (cos x)): sin/cos kernels %n, the bare div per door r
  static float32_t _rs_tan(float32_t x) {
    float32_t s = _rs_sin(x), c = _rs_cos(x);
    softfloat_roundingMode = _math_rnd;
    float32_t r = f32_div(s, c);
    softfloat_roundingMode = softfloat_round_near_even;
    return r;
  }

/* @rs sqt -- math.hoon ++rs ++sqt = (sqt:^rs x): correctly-rounded f32 sqrt. */
  static float32_t _rs_sqt(float32_t x) {
    union sing r; r.s = f32_sqrt(x);
    if ( !f32_eq(r.s, r.s) ) r.c = _RS_QNAN;        // _nan_unify
    return r.s;
  }

/* @rs log/log-2/log-10 -- math.hoon ++rs ++log/++lr/++log-2/++log-10
**   x = 2^e * m, m in [sqrt(1/2),sqrt(2)); log(1+f) via atanh series (deg-4).
*/
  //  +lr: finite positive x -> *ef = e as @rs, *l1 = log(mantissa)
  static void _rs_lr(float32_t x, float32_t* ef, float32_t* l1) {
    static const uint32_t cs[5] = {
      0x3eaaaaabU, 0x3e4ccccdU, 0x3e124925U, 0x3de38e39U, 0x3dba2e8cU };
    union sing xb, m, f, s, z, p2, r, ll, efa, c, one, half;
    one.c = 0x3f800000U; half.c = 0x3f000000U;
    xb.s = x;
    int sub = (((xb.c >> 23) & 0xffU) == 0);
    if ( sub ) xb.s = f32_mul(x, _rs_bits(0x4b800000U).s);     // *2^24
    int32_t ae = sub ? -24 : 0;
    int32_t e = (int32_t)((xb.c >> 23) & 0xffU) - 127;
    m.c = (xb.c & 0x7fffffU) | 0x3f800000U;
    if ( !f32_lt(m.s, _rs_bits(0x3fb504f3U).s) ) {            // m >= sqrt(2)
      m.s = f32_mul(m.s, half.s); e += 1;
    }
    e += ae;
    f.s = f32_sub(m.s, one.s);
    s.s = f32_div(f.s, f32_add(m.s, one.s));
    z.s = f32_mul(s.s, s.s);
    p2.c = 0; for ( int i = 5; i-- != 0; ) { c.c = cs[i]; p2.s = f32_add(f32_mul(p2.s, z.s), c.s); }
    r.s = f32_mul(f32_add(z.s, z.s), p2.s);
    ll.s = f32_sub(f.s, f32_mul(s.s, f32_sub(f.s, r.s)));
    efa.s = ui32_to_f32((uint32_t)(e < 0 ? -e : e));
    *ef = (e >= 0) ? efa.s : _rs_neg(efa.s);
    *l1 = ll.s;
  }
  //  shared guards for log/log-2/log-10; returns 1 (and sets *g) on a special case
  static int _rs_log_guard(float32_t x, float32_t* g) {
    union sing r0; r0.s = x;
    if ( !f32_eq(x, x) )    { r0.c = _RS_QNAN; *g = r0.s; return 1; }       // NaN
    if ( r0.c == _RS_PINF ) { *g = x; return 1; }                          // +inf -> inf
    if ( (r0.c == 0)||(r0.c == 0x80000000U) ) { r0.c = _RS_NINF; *g = r0.s; return 1; }  // +-0 -> -inf
    if ( (r0.c >> 31) == 1 ){ r0.c = _RS_QNAN; *g = r0.s; return 1; }       // x<0 -> NaN
    return 0;
  }
  static float32_t _rs_log(float32_t x) {
    union sing g, ef, l1, hi, lo;
    if ( _rs_log_guard(x, &g.s) ) return g.s;
    _rs_lr(x, &ef.s, &l1.s);
    hi.s = f32_mul(ef.s, _rs_bits(0x3f317200U).s);                // e*ln2hi
    lo.s = f32_mul(ef.s, _rs_bits(0x35bfbe8eU).s);                // e*ln2lo
    return f32_add(hi.s, f32_add(l1.s, lo.s));
  }
  static float32_t _rs_log2(float32_t x) {
    union sing g, ef, l1;
    if ( _rs_log_guard(x, &g.s) ) return g.s;
    _rs_lr(x, &ef.s, &l1.s);
    return f32_add(ef.s, f32_mul(l1.s, _rs_bits(0x3fb8aa3bU).s)); // e + lm/ln2
  }
  static float32_t _rs_log10(float32_t x) {
    union sing g, ef, l1;
    if ( _rs_log_guard(x, &g.s) ) return g.s;
    _rs_lr(x, &ef.s, &l1.s);
    return f32_add(f32_mul(ef.s, _rs_bits(0x3e9a209bU).s),        // e*log10(2)
                   f32_mul(l1.s, _rs_bits(0x3ede5bd9U).s));       // + lm/ln10
  }

/* @rs cbt -- math.hoon ++rs ++cbt = sign(x) * exp(log|x| / 3). */
  static float32_t _rs_cbt(float32_t x) {
    union sing r0, ax, r;
    r0.s = x;
    if ( !f32_eq(x, x) )                       { return x; }                      // NaN
    if ( (r0.c == 0)||(r0.c == 0x80000000U) )  { return x; }                      // +-0
    ax.c = r0.c & 0x7fffffffU;
    r.s = _rs_exp(f32_mul(_rs_log(ax.s), _rs_bits(0x3eaaaaabU).s));               // exp(log|x|/3)
    return ((r0.c >> 31) == 1) ? _rs_neg(r.s) : r.s;
  }

/* @rs asin/acos -- math.hoon ++rs ++asin/++acos/++rs-ainv
**   rational kernel R(t) = t*P2(t)/(1 + c*t); sqt = _rs_sqt.
*/
  static float32_t _rs_ainv_rr(float32_t t) {
    static const uint32_t ps[3] = { 0x3e2aaa75U, 0xbd2f13baU, 0xbc0dd36bU };
    union sing pp, c, one;
    one.c = 0x3f800000U;
    pp.c = 0; for ( int i = 3; i-- != 0; ) { c.c = ps[i]; pp.s = f32_add(f32_mul(pp.s, t), c.s); }
    return f32_div(f32_mul(t, pp.s),
                   f32_add(one.s, f32_mul(t, _rs_bits(0xbf34e5aeU).s)));
  }
  static float32_t _rs_asin(float32_t x) {
    union sing r0, ax, t, w, r, s, res, half, one, two, pio2h, pio2l, pio4;
    half.c=0x3f000000U; one.c=0x3f800000U; two.c=0x40000000U;
    pio2h.c=0x3fc90fdbU; pio2l.c=0xb33bbd2eU; pio4.c=0x3f490fdbU;
    r0.s = x;
    if ( !f32_eq(x, x) )       { r0.c = _RS_QNAN; return r0.s; }   // NaN
    uint32_t sgn = r0.c >> 31;
    ax.c = r0.c & 0x7fffffffU;
    if ( f32_lt(one.s, ax.s) ) { r0.c = _RS_QNAN; return r0.s; }   // |x|>1 -> NaN
    if ( ax.c == one.c )                                          // |x|==1
      return f32_add(f32_mul(x, pio2h.s), f32_mul(x, pio2l.s));
    if ( f32_lt(ax.s, half.s) ) {                                // |x|<0.5
      if ( f32_lt(ax.s, _rs_bits(0x39800000U).s) ) return x;      // tiny
      t.s = f32_mul(x, x);
      return f32_add(x, f32_mul(x, _rs_ainv_rr(t.s)));
    }
    w.s = f32_sub(one.s, ax.s);
    t.s = f32_mul(w.s, half.s);
    r.s = _rs_ainv_rr(t.s);
    s.s = _rs_sqt(t.s);
    if ( f32_le(_rs_bits(0x3f79999aU).s, ax.s) ) {                // near 1
      res.s = f32_sub(pio2h.s, f32_sub(f32_mul(two.s, f32_add(s.s, f32_mul(s.s, r.s))), pio2l.s));
      return (sgn == 1) ? _rs_neg(res.s) : res.s;
    }
    { union sing df, cc, p2, q2;
      df.c = s.c & 0xfffff000U;
      cc.s = f32_div(f32_sub(t.s, f32_mul(df.s, df.s)), f32_add(s.s, df.s));
      p2.s = f32_sub(f32_mul(two.s, f32_mul(s.s, r.s)), f32_sub(pio2l.s, f32_mul(two.s, cc.s)));
      q2.s = f32_sub(pio4.s, f32_mul(two.s, df.s));
      res.s = f32_sub(pio4.s, f32_sub(p2.s, q2.s));
      return (sgn == 1) ? _rs_neg(res.s) : res.s;
    }
  }
  static float32_t _rs_acos(float32_t x) {
    union sing r0, ax, z, s, r, w, half, one, two, pi, pio2h, pio2l;
    half.c=0x3f000000U; one.c=0x3f800000U; two.c=0x40000000U;
    pi.c=0x40490fdbU; pio2h.c=0x3fc90fdbU; pio2l.c=0xb33bbd2eU;
    r0.s = x;
    if ( !f32_eq(x, x) )       { r0.c = _RS_QNAN; return r0.s; }   // NaN
    uint32_t neg = r0.c >> 31;
    ax.c = r0.c & 0x7fffffffU;
    if ( f32_lt(one.s, ax.s) ) { r0.c = _RS_QNAN; return r0.s; }   // |x|>1 -> NaN
    if ( ax.c == one.c ) {                                        // |x|==1
      if ( neg == 0 ) { union sing z0; z0.c = 0; return z0.s; }    // 1 -> 0
      return f32_add(pi.s, f32_mul(two.s, pio2l.s));               // -1 -> pi
    }
    if ( f32_lt(ax.s, half.s) ) {                                // |x|<0.5
      if ( f32_lt(ax.s, _rs_bits(0x32800000U).s) ) return pio2h.s; // tiny -> pi/2
      z.s = f32_mul(x, x);
      r.s = _rs_ainv_rr(z.s);
      return f32_sub(pio2h.s, f32_sub(x, f32_sub(pio2l.s, f32_mul(x, r.s))));
    }
    if ( neg == 1 ) {                                            // x <= -0.5
      z.s = f32_mul(f32_add(one.s, x), half.s);
      s.s = _rs_sqt(z.s);
      r.s = _rs_ainv_rr(z.s);
      w.s = f32_sub(f32_mul(r.s, s.s), pio2l.s);
      return f32_sub(pi.s, f32_mul(two.s, f32_add(s.s, w.s)));
    }
    z.s = f32_mul(f32_sub(one.s, x), half.s);                    // x >= 0.5
    s.s = _rs_sqt(z.s);
    r.s = _rs_ainv_rr(z.s);
    return f32_mul(two.s, f32_add(s.s, f32_mul(s.s, r.s)));
  }

/* @rs atan/atan2 -- math.hoon ++rs ++atan/++rs-atan/++atan2
**   fdlibm breakpoint reduction (7/16,11/16,19/16,39/16) + deg-4 minimax.
*/
  static float32_t _rs_atan(float32_t x) {
    static const uint32_t at[5] = {
      0x3eaaaaa9U, 0xbe4cca98U, 0x3e11f50dU, 0xbdda1247U, 0x3d7cac25U };
    union sing r0, ax, xr, hi, lo, z, sp, s, res, one, two, ohf, c;
    r0.s = x;
    if ( !f32_eq(x, x) )    { r0.c = _RS_QNAN; return r0.s; }       // NaN
    if ( r0.c == _RS_PINF ) { r0.c = 0x3fc90fdbU; return r0.s; }    // +inf -> pi/2
    if ( r0.c == _RS_NINF ) { r0.c = 0xbfc90fdbU; return r0.s; }    // -inf -> -pi/2
    if ( (r0.c == 0)||(r0.c == 0x80000000U) ) { return x; }         // +-0
    uint32_t neg = r0.c >> 31;
    ax.c = r0.c & 0x7fffffffU;
    one.c=0x3f800000U; two.c=0x40000000U; ohf.c=0x3fc00000U;
    int dir = 0;
    if ( f32_lt(ax.s, _rs_bits(0x3ee00000U).s) ) {                 // |x| < 7/16
      xr.s = ax.s; hi.c = 0; lo.c = 0; dir = 1;
    } else if ( f32_lt(ax.s, _rs_bits(0x3f300000U).s) ) {          // < 11/16
      xr.s = f32_div(f32_sub(f32_add(ax.s, ax.s), one.s), f32_add(two.s, ax.s));
      hi.c = 0x3eed6338U; lo.c = 0x31ac376aU;                      // atan(0.5)
    } else if ( f32_lt(ax.s, _rs_bits(0x3f980000U).s) ) {          // < 19/16
      xr.s = f32_div(f32_sub(ax.s, one.s), f32_add(ax.s, one.s));
      hi.c = 0x3f490fdbU; lo.c = 0xb2bbbd2eU;                      // pi/4
    } else if ( f32_lt(ax.s, _rs_bits(0x401c0000U).s) ) {          // < 39/16
      xr.s = f32_div(f32_sub(ax.s, ohf.s), f32_add(one.s, f32_mul(ohf.s, ax.s)));
      hi.c = 0x3f7b985fU; lo.c = 0xb2d7e096U;                      // atan(1.5)
    } else {                                                       // -1/x
      xr.s = f32_div(_rs_bits(0xbf800000U).s, ax.s);
      hi.c = 0x3fc90fdbU; lo.c = 0xb33bbd2eU;                      // pi/2
    }
    z.s = f32_mul(xr.s, xr.s);
    sp.c = 0; for ( int i = 5; i-- != 0; ) { c.c = at[i]; sp.s = f32_add(f32_mul(sp.s, z.s), c.s); }
    s.s = f32_mul(z.s, sp.s);
    if ( dir ) res.s = f32_sub(xr.s, f32_mul(xr.s, s.s));
    else       res.s = f32_sub(hi.s, f32_sub(f32_sub(f32_mul(xr.s, s.s), lo.s), xr.s));
    return (neg == 1) ? _rs_neg(res.s) : res.s;
  }
  //  bare door ops (div/add/sub/mul) round per _math_rnd; atan kernel is %n.
  static float32_t _rs_atan2(float32_t y, float32_t x) {
    union sing xb, pi, two, zero, mone, q, a, r;
    zero.c = 0; pi.c = 0x40490fdbU; two.c = 0x40000000U; mone.c = 0xbf800000U;
    xb.s = x;
    if ( f32_lt(zero.s, x) ) {                                     // x>0: atan(div y x)
      softfloat_roundingMode = _math_rnd; q.s = f32_div(y, x);
      softfloat_roundingMode = softfloat_round_near_even; return _rs_atan(q.s);
    }
    if ( f32_lt(x, zero.s) && f32_le(zero.s, y) ) {                // x<0,y>=0: add(atan,pi)
      softfloat_roundingMode = _math_rnd; q.s = f32_div(y, x);
      softfloat_roundingMode = softfloat_round_near_even; a.s = _rs_atan(q.s);
      softfloat_roundingMode = _math_rnd; r.s = f32_add(a.s, pi.s);
      softfloat_roundingMode = softfloat_round_near_even; return r.s;
    }
    if ( f32_lt(x, zero.s) && f32_lt(y, zero.s) ) {                // x<0,y<0: sub(atan,pi)
      softfloat_roundingMode = _math_rnd; q.s = f32_div(y, x);
      softfloat_roundingMode = softfloat_round_near_even; a.s = _rs_atan(q.s);
      softfloat_roundingMode = _math_rnd; r.s = f32_sub(a.s, pi.s);
      softfloat_roundingMode = softfloat_round_near_even; return r.s;
    }
    if ( (xb.c == 0) && f32_lt(zero.s, y) ) {                      // x==+0,y>0: div(pi,2)
      softfloat_roundingMode = _math_rnd; r.s = f32_div(pi.s, two.s);
      softfloat_roundingMode = softfloat_round_near_even; return r.s;
    }
    if ( (xb.c == 0) && f32_lt(y, zero.s) ) {                      // x==+0,y<0: mul(-1,div(pi,2))
      softfloat_roundingMode = _math_rnd;
      r.s = f32_mul(mone.s, f32_div(pi.s, two.s));
      softfloat_roundingMode = softfloat_round_near_even; return r.s;
    }
    return zero.s;
  }

/* @rs pow/pow-n -- math.hoon ++rs ++pow/++pow-n */
  static float32_t _rs_pow_n(float32_t x, float32_t n) {
    union sing nn, p, one, two;
    one.c = 0x3f800000U; two.c = 0x40000000U;
    nn.s = n;
    if ( nn.c == 0 ) return one.s;                 // n == +0 -> 1
    softfloat_roundingMode = _math_rnd;            // bare mul/sub round per door r
    p.s = x;
    while ( !f32_lt(n, two.s) ) { p.s = f32_mul(p.s, x); n = f32_sub(n, one.s); }
    softfloat_roundingMode = softfloat_round_near_even;
    return p.s;
  }
  static float32_t _rs_pow(float32_t x, float32_t n) {
    union sing nn, ni, zero, lg, prod;
    zero.c = 0; nn.s = n;
    ni.s = i32_to_f32(f32_to_i32(n, softfloat_round_near_even, 0));   // san (need (toi n)) (mode-indep)
    if ( (nn.c == ni.c) && f32_lt(zero.s, n) )                        // positive integer
      return _rs_pow_n(x, ni.s);
    lg.s = _rs_log(x);                                               // %n kernel
    softfloat_roundingMode = _math_rnd;                             // bare mul per door r
    prod.s = f32_mul(n, lg.s);
    softfloat_roundingMode = softfloat_round_near_even;
    return _rs_exp(prod.s);                                          // %n kernel: exp(n*log x)
  }

/* ===================================================================
** @rh (half-precision) cores -- math.hoon ++rh.  Every arm is
**   narrow-sh(rs_fn(widen-hs(x))): widen f16->f32 (exact = f16_to_f32),
**   compute in @rs, narrow f32->f16 (= f32_to_f16) honoring the door's r.
** The rs call inherits r too (composites honor it; kernels stay %n).  No new
** algorithm cores -- the @rs cores are reused.  Verified: the Hoon narrow-sh
** is bit-exact to f32_to_f16 in all four rounding modes.
** =================================================================== */

  union half {
    float16_t h;
    uint16_t  c;
  };

  //  widen f16->f32 (exact), run an @rs core, narrow f32->f16 per _math_rnd.
  //  _math_rnd is the door's r (set by the wrapper): the rs composites honor it
  //  and the narrow rounds by it; the rs kernels run at near-even.
  static inline float16_t _rh_1(float16_t x, float32_t (*fun)(float32_t)) {
    softfloat_roundingMode = softfloat_round_near_even;
    float32_t v = fun(f16_to_f32(x));
    softfloat_roundingMode = _math_rnd;
    return f32_to_f16(v);
  }
  static inline float16_t _rh_2(float16_t x, float16_t n,
                                float32_t (*fun)(float32_t, float32_t)) {
    softfloat_roundingMode = softfloat_round_near_even;
    float32_t v = fun(f16_to_f32(x), f16_to_f32(n));
    softfloat_roundingMode = _math_rnd;
    return f32_to_f16(v);
  }
  static float16_t _rh_exp(float16_t x)   { return _rh_1(x, _rs_exp); }
  static float16_t _rh_log(float16_t x)   { return _rh_1(x, _rs_log); }
  static float16_t _rh_sin(float16_t x)   { return _rh_1(x, _rs_sin); }
  static float16_t _rh_cos(float16_t x)   { return _rh_1(x, _rs_cos); }
  static float16_t _rh_tan(float16_t x)   { return _rh_1(x, _rs_tan); }
  static float16_t _rh_atan(float16_t x)  { return _rh_1(x, _rs_atan); }
  static float16_t _rh_asin(float16_t x)  { return _rh_1(x, _rs_asin); }
  static float16_t _rh_acos(float16_t x)  { return _rh_1(x, _rs_acos); }
  static float16_t _rh_sqt(float16_t x)   { return _rh_1(x, _rs_sqt); }
  static float16_t _rh_cbt(float16_t x)   { return _rh_1(x, _rs_cbt); }
  static float16_t _rh_log2(float16_t x)  { return _rh_1(x, _rs_log2); }
  static float16_t _rh_log10(float16_t x) { return _rh_1(x, _rs_log10); }
  static float16_t _rh_atan2(float16_t y, float16_t x) { return _rh_2(y, x, _rs_atan2); }
  static float16_t _rh_pow(float16_t x, float16_t n)   { return _rh_2(x, n, _rs_pow); }
  static float16_t _rh_pow_n(float16_t x, float16_t n) { return _rh_2(x, n, _rs_pow_n); }

/* ===================================================================
** @rq (quad, 128-bit) cores -- math.hoon ++rq.  Native f128 algorithms (the
** widest type, no delegation): same reductions as @rd, higher-degree minimax
** in float128_t.  Marshalling reads TWO chubs into float128_t.v[0..1], so it is
** word-size-agnostic -- the chub ABI sidesteps the c3_w*[n=2|4] divergence of
** the old rq.c.  Composite arms honor the door's r via _math_rnd.
** =================================================================== */

  union quad {
    float128_t q;
    c3_d       w[2];        // w[0] = v[0] = low 64 bits, w[1] = v[1] = high 64
  };
  static inline float128_t _rq_bits(c3_d hi, c3_d lo) {
    union quad u; u.w[0] = lo; u.w[1] = hi; return u.q;
  }
  //  by-value wrappers over the pointer-based f128M_* ops (this SoftFloat build
  //  has no by-value f128_*).  Compiler inlines these; keeps the cores readable.
  static inline float128_t _rqm(float128_t a, float128_t b) { float128_t r; f128M_mul(&a,&b,&r); return r; }
  static inline float128_t _rqa(float128_t a, float128_t b) { float128_t r; f128M_add(&a,&b,&r); return r; }
  static inline float128_t _rqs(float128_t a, float128_t b) { float128_t r; f128M_sub(&a,&b,&r); return r; }
  static inline float128_t _rqd(float128_t a, float128_t b) { float128_t r; f128M_div(&a,&b,&r); return r; }
  static inline float128_t _rqq(float128_t a)               { float128_t r; f128M_sqrt(&a,&r);  return r; }
  static inline int _rqeq(float128_t a, float128_t b) { return f128M_eq(&a,&b); }
  static inline int _rqlt(float128_t a, float128_t b) { return f128M_lt(&a,&b); }
  static inline int _rqle(float128_t a, float128_t b) { return f128M_le(&a,&b); }
  static inline c3_ds _rqtoi(float128_t a, int m) { return (c3_ds)f128M_to_i64(&a, (uint_fast8_t)m, 0); }
  static inline float128_t _rqi64(c3_ds n) { float128_t r; i64_to_f128M(n, &r); return r; }
  static const c3_d _RQ_QNAN_HI = 0x7fff800000000000ULL;
  static const c3_d _RQ_PINF_HI = 0x7fff000000000000ULL;
  static const c3_d _RQ_NINF_HI = 0xffff000000000000ULL;
  static inline float128_t _rq_neg(float128_t a) {        // (sub .0 a)
    return _rqs(_rq_bits(0,0), a);
  }

/* @rq exp -- math.hoon ++rq ++exp
**   Cody-Waite x=k*ln2+r, exp=2^k*P(r), P a degree-24 minimax (f128).
*/
  //  pow2(j) = 2^j as f128 bits: exp field (bits 112-126) = j+16383
  static inline float128_t _rq_pow2(c3_ds j) {
    union quad u; u.w[0] = 0; u.w[1] = ((c3_d)(j + 16383)) << 48; return u.q;
  }
  static inline float128_t _rq_scale2(float128_t p, c3_ds k) {
    if ( (k - 16384) >= 0 )
      return _rqm(_rqm(p, _rq_pow2(16383)), _rq_pow2(k - 16383));
    if ( !((k + 16382) >= 0) )
      return _rqm(_rqm(p, _rq_pow2(k + 112)), _rq_pow2(-112));
    return _rqm(p, _rq_pow2(k));
  }
  static float128_t _rq_exp(float128_t x) {
    //  degree-24 minimax coeffs c0..c24 {lo, hi} (math.hoon ++rq ++exp)
    static const c3_d cs[25][2] = {
      {0x0000000000000000ULL,0x3fff000000000000ULL},{0x0000000000000000ULL,0x3fff000000000000ULL},
      {0x0000000000000000ULL,0x3ffe000000000000ULL},{0x5555555555555555ULL,0x3ffc555555555555ULL},
      {0x5555555555555555ULL,0x3ffa555555555555ULL},{0x1111111111111111ULL,0x3ff8111111111111ULL},
      {0x6c16c16c16c16c17ULL,0x3ff56c16c16c16c1ULL},{0xa01a01a01a01a3e8ULL,0x3ff2a01a01a01a01ULL},
      {0xa01a01a01a01a146ULL,0x3fefa01a01a01a01ULL},{0x38faac1c88a5a526ULL,0x3fec71de3a556c73ULL},
      {0xc72ef016d3d6e867ULL,0x3fe927e4fb7789f5ULL},{0x38fe748363c46e8bULL,0x3fe5ae64567f544eULL},
      {0x7b544dab18f475c5ULL,0x3fe21eed8eff8d89ULL},{0x97c9f3aebabb2423ULL,0x3fde6124613a86d0ULL},
      {0xd20b83c7f94d17d8ULL,0x3fda93974a8c07c9ULL},{0xf5f4284f0d74f9e7ULL,0x3fd6ae7f3e733b81ULL},
      {0xf417b4d27c5f92a9ULL,0x3fd2ae7f3e733b81ULL},{0x6a419e674779c97cULL,0x3fce952c77030a99ULL},
      {0x0466ff8c8b42b3dfULL,0x3fca6827863b97b5ULL},{0x874b7a686d819241ULL,0x3fc62f49b469f892ULL},
      {0xbb3b32a11bb5f139ULL,0x3fc1e542ba427463ULL},{0xc93890ff9ab55cbbULL,0x3fbd71b8db9f7f73ULL},
      {0x6efc0717eae785a1ULL,0x3fb90ce38aab7bd7ULL},{0xcb3f4f7edfaa2666ULL,0x3fb47693274bab2aULL},
      {0x61cb0e23655d47cbULL,0x3faff3629154e0a7ULL},
    };
    union quad r0; r0.q = x;
    if ( !_rqeq(x, x) )                       return _rq_bits(_RQ_QNAN_HI, 0);   // NaN
    if ( r0.w[1]==_RQ_PINF_HI && r0.w[0]==0 ) return x;                          // +inf
    if ( r0.w[1]==_RQ_NINF_HI && r0.w[0]==0 ) return _rq_bits(0,0);              // -inf -> 0

    float128_t log2e = _rq_bits(0x3fff71547652b82fULL, 0xe1777d0ffda0d23aULL);
    float128_t ln2hi = _rq_bits(0x3ffe62e42fefa39eULL, 0xf35793c800000000ULL);
    float128_t ln2lo = _rq_bits(0xbfad319ff0342542ULL, 0xfc32f366359d274aULL);

    c3_ds k = _rqtoi(_rqm(x, log2e), softfloat_round_near_even);
    if ( (k - 16385) >= 0 )    return _rq_bits(_RQ_PINF_HI, 0);                  // overflow -> inf
    if ( !((k + 16494) >= 0) ) return _rq_bits(0, 0);                           // underflow -> 0

    float128_t ka = _rqi64((c3_ds)(k < 0 ? -k : k));
    float128_t kf = (k >= 0) ? ka : _rq_neg(ka);
    float128_t rr = _rqs( _rqs(x, _rqm(kf, ln2hi)), _rqm(kf, ln2lo) );

    float128_t p = _rq_bits(0,0);
    for ( int i = 25; i-- != 0; )          // Horner over flop(cs): c24..c0
      p = _rqa(_rqm(p, rr), _rq_bits(cs[i][1], cs[i][0]));
    return _rq_scale2(p, k);
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
      _math_rnd = _rnd_of(u3r_at(60, cor));      // door rounding r (for composite arms)
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
    if ( c3n == u3r_mean(cor, u3x_sam_2, &y, u3x_sam_3, &x, 0) ||
         c3n == u3ud(y) || c3n == u3ud(x) ) {
      return u3m_bail(c3__exit);
    }
    _math_rnd = _rnd_of(u3r_at(60, cor));         // door rounding r
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
    if ( c3n == u3r_mean(cor, u3x_sam_2, &x, u3x_sam_3, &n, 0) ||
         c3n == u3ud(x) || c3n == u3ud(n) ) {
      return u3m_bail(c3__exit);
    }
    { union doub xx, nn, e;
      softfloat_roundingMode = softfloat_round_near_even;
      _math_rnd = _rnd_of(u3r_at(60, cor));      // door rounding r
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

/* @rs ABI wrappers.  @rs is a 32-bit atom: read the low 32 bits of the chub,
** write the 32-bit result as a chub (high bits zero -> normalizes to a 32-bit
** atom).  Chub I/O keeps this word-size-agnostic, like the @rd wrappers.
*/
  static inline float32_t _rs_in(u3_atom a) {
    union sing s; s.c = (uint32_t)u3r_chub(0, a); return s.s;
  }
  static inline u3_noun _rs_out(float32_t v) {
    union sing s; s.s = v; { c3_d out = (c3_d)s.c; return u3i_chubs(1, &out); }
  }
  static u3_noun _rs_jet(u3_noun cor, float32_t (*fun)(float32_t)) {
    u3_noun x = u3r_at(u3x_sam, cor);
    if ( u3_none == x || c3n == u3ud(x) ) return u3m_bail(c3__exit);
    softfloat_roundingMode = softfloat_round_near_even;
    _math_rnd = _rnd_of(u3r_at(60, cor));        // door rounding r (for @rs tan)
    return _rs_out(fun(_rs_in(x)));
  }
  static u3_noun _rs_jet2(u3_noun cor, float32_t (*fun)(float32_t, float32_t)) {
    u3_noun x, n;
    if ( c3n == u3r_mean(cor, u3x_sam_2, &x, u3x_sam_3, &n, 0) ||
         c3n == u3ud(x) || c3n == u3ud(n) ) {
      return u3m_bail(c3__exit);
    }
    softfloat_roundingMode = softfloat_round_near_even;
    _math_rnd = _rnd_of(u3r_at(60, cor));        // door rounding r
    return _rs_out(fun(_rs_in(x), _rs_in(n)));
  }

  u3_noun u3qi_rs_exp(u3_atom a)   { softfloat_roundingMode=softfloat_round_near_even; return _rs_out(_rs_exp(_rs_in(a))); }
  u3_noun u3wi_rs_exp(u3_noun cor) { return _rs_jet(cor, _rs_exp); }
  u3_noun u3qi_rs_log(u3_atom a)   { softfloat_roundingMode=softfloat_round_near_even; return _rs_out(_rs_log(_rs_in(a))); }
  u3_noun u3wi_rs_log(u3_noun cor) { return _rs_jet(cor, _rs_log); }
  u3_noun u3qi_rs_sin(u3_atom a)   { softfloat_roundingMode=softfloat_round_near_even; return _rs_out(_rs_sin(_rs_in(a))); }
  u3_noun u3wi_rs_sin(u3_noun cor) { return _rs_jet(cor, _rs_sin); }
  u3_noun u3qi_rs_cos(u3_atom a)   { softfloat_roundingMode=softfloat_round_near_even; return _rs_out(_rs_cos(_rs_in(a))); }
  u3_noun u3wi_rs_cos(u3_noun cor) { return _rs_jet(cor, _rs_cos); }
  u3_noun u3qi_rs_tan(u3_atom a)   { softfloat_roundingMode=softfloat_round_near_even; return _rs_out(_rs_tan(_rs_in(a))); }
  u3_noun u3wi_rs_tan(u3_noun cor) { return _rs_jet(cor, _rs_tan); }
  u3_noun u3qi_rs_atan(u3_atom a)  { softfloat_roundingMode=softfloat_round_near_even; return _rs_out(_rs_atan(_rs_in(a))); }
  u3_noun u3wi_rs_atan(u3_noun cor){ return _rs_jet(cor, _rs_atan); }
  u3_noun u3qi_rs_asin(u3_atom a)  { softfloat_roundingMode=softfloat_round_near_even; return _rs_out(_rs_asin(_rs_in(a))); }
  u3_noun u3wi_rs_asin(u3_noun cor){ return _rs_jet(cor, _rs_asin); }
  u3_noun u3qi_rs_acos(u3_atom a)  { softfloat_roundingMode=softfloat_round_near_even; return _rs_out(_rs_acos(_rs_in(a))); }
  u3_noun u3wi_rs_acos(u3_noun cor){ return _rs_jet(cor, _rs_acos); }
  u3_noun u3qi_rs_sqt(u3_atom a)   { softfloat_roundingMode=softfloat_round_near_even; return _rs_out(_rs_sqt(_rs_in(a))); }
  u3_noun u3wi_rs_sqt(u3_noun cor) { return _rs_jet(cor, _rs_sqt); }
  u3_noun u3qi_rs_cbt(u3_atom a)   { softfloat_roundingMode=softfloat_round_near_even; return _rs_out(_rs_cbt(_rs_in(a))); }
  u3_noun u3wi_rs_cbt(u3_noun cor) { return _rs_jet(cor, _rs_cbt); }
  u3_noun u3qi_rs_log2(u3_atom a)  { softfloat_roundingMode=softfloat_round_near_even; return _rs_out(_rs_log2(_rs_in(a))); }
  u3_noun u3wi_rs_log2(u3_noun cor){ return _rs_jet(cor, _rs_log2); }
  u3_noun u3qi_rs_log10(u3_atom a) { softfloat_roundingMode=softfloat_round_near_even; return _rs_out(_rs_log10(_rs_in(a))); }
  u3_noun u3wi_rs_log10(u3_noun cor){ return _rs_jet(cor, _rs_log10); }

  u3_noun u3qi_rs_atan2(u3_atom y, u3_atom x) { softfloat_roundingMode=softfloat_round_near_even; return _rs_out(_rs_atan2(_rs_in(y), _rs_in(x))); }
  u3_noun u3wi_rs_atan2(u3_noun cor){ return _rs_jet2(cor, _rs_atan2); }
  u3_noun u3qi_rs_pow(u3_atom x, u3_atom n)   { softfloat_roundingMode=softfloat_round_near_even; return _rs_out(_rs_pow(_rs_in(x), _rs_in(n))); }
  u3_noun u3wi_rs_pow(u3_noun cor) { return _rs_jet2(cor, _rs_pow); }
  u3_noun u3qi_rs_pow_n(u3_atom x, u3_atom n) { softfloat_roundingMode=softfloat_round_near_even; return _rs_out(_rs_pow_n(_rs_in(x), _rs_in(n))); }
  u3_noun u3wi_rs_pow_n(u3_noun cor){ return _rs_jet2(cor, _rs_pow_n); }

/* @rh ABI wrappers.  @rh is a 16-bit atom: read the low 16 bits of the chub,
** write the 16-bit result via chub (high bits zero -> normalizes).  Same
** word-agnostic chub I/O as @rd/@rs.  Wrappers set _math_rnd from the door's r
** (axis 60); the cores apply it to the rs composite ops and the f16 narrow.
*/
  static inline float16_t _rh_in(u3_atom a) {
    union half s; s.c = (uint16_t)u3r_chub(0, a); return s.h;
  }
  static inline u3_noun _rh_out(float16_t v) {
    union half s; s.h = v; { c3_d out = (c3_d)s.c; return u3i_chubs(1, &out); }
  }
  static u3_noun _rh_jet(u3_noun cor, float16_t (*fun)(float16_t)) {
    u3_noun x = u3r_at(u3x_sam, cor);
    if ( u3_none == x || c3n == u3ud(x) ) return u3m_bail(c3__exit);
    _math_rnd = _rnd_of(u3r_at(60, cor));
    return _rh_out(fun(_rh_in(x)));
  }
  static u3_noun _rh_jet2(u3_noun cor, float16_t (*fun)(float16_t, float16_t)) {
    u3_noun x, n;
    if ( c3n == u3r_mean(cor, u3x_sam_2, &x, u3x_sam_3, &n, 0) ||
         c3n == u3ud(x) || c3n == u3ud(n) ) {
      return u3m_bail(c3__exit);
    }
    _math_rnd = _rnd_of(u3r_at(60, cor));
    return _rh_out(fun(_rh_in(x), _rh_in(n)));
  }

  u3_noun u3qi_rh_exp(u3_atom a)   { _math_rnd=softfloat_round_near_even; return _rh_out(_rh_exp(_rh_in(a))); }
  u3_noun u3wi_rh_exp(u3_noun cor) { return _rh_jet(cor, _rh_exp); }
  u3_noun u3qi_rh_log(u3_atom a)   { _math_rnd=softfloat_round_near_even; return _rh_out(_rh_log(_rh_in(a))); }
  u3_noun u3wi_rh_log(u3_noun cor) { return _rh_jet(cor, _rh_log); }
  u3_noun u3qi_rh_sin(u3_atom a)   { _math_rnd=softfloat_round_near_even; return _rh_out(_rh_sin(_rh_in(a))); }
  u3_noun u3wi_rh_sin(u3_noun cor) { return _rh_jet(cor, _rh_sin); }
  u3_noun u3qi_rh_cos(u3_atom a)   { _math_rnd=softfloat_round_near_even; return _rh_out(_rh_cos(_rh_in(a))); }
  u3_noun u3wi_rh_cos(u3_noun cor) { return _rh_jet(cor, _rh_cos); }
  u3_noun u3qi_rh_tan(u3_atom a)   { _math_rnd=softfloat_round_near_even; return _rh_out(_rh_tan(_rh_in(a))); }
  u3_noun u3wi_rh_tan(u3_noun cor) { return _rh_jet(cor, _rh_tan); }
  u3_noun u3qi_rh_atan(u3_atom a)  { _math_rnd=softfloat_round_near_even; return _rh_out(_rh_atan(_rh_in(a))); }
  u3_noun u3wi_rh_atan(u3_noun cor){ return _rh_jet(cor, _rh_atan); }
  u3_noun u3qi_rh_asin(u3_atom a)  { _math_rnd=softfloat_round_near_even; return _rh_out(_rh_asin(_rh_in(a))); }
  u3_noun u3wi_rh_asin(u3_noun cor){ return _rh_jet(cor, _rh_asin); }
  u3_noun u3qi_rh_acos(u3_atom a)  { _math_rnd=softfloat_round_near_even; return _rh_out(_rh_acos(_rh_in(a))); }
  u3_noun u3wi_rh_acos(u3_noun cor){ return _rh_jet(cor, _rh_acos); }
  u3_noun u3qi_rh_sqt(u3_atom a)   { _math_rnd=softfloat_round_near_even; return _rh_out(_rh_sqt(_rh_in(a))); }
  u3_noun u3wi_rh_sqt(u3_noun cor) { return _rh_jet(cor, _rh_sqt); }
  u3_noun u3qi_rh_cbt(u3_atom a)   { _math_rnd=softfloat_round_near_even; return _rh_out(_rh_cbt(_rh_in(a))); }
  u3_noun u3wi_rh_cbt(u3_noun cor) { return _rh_jet(cor, _rh_cbt); }
  u3_noun u3qi_rh_log2(u3_atom a)  { _math_rnd=softfloat_round_near_even; return _rh_out(_rh_log2(_rh_in(a))); }
  u3_noun u3wi_rh_log2(u3_noun cor){ return _rh_jet(cor, _rh_log2); }
  u3_noun u3qi_rh_log10(u3_atom a) { _math_rnd=softfloat_round_near_even; return _rh_out(_rh_log10(_rh_in(a))); }
  u3_noun u3wi_rh_log10(u3_noun cor){ return _rh_jet(cor, _rh_log10); }

  u3_noun u3qi_rh_atan2(u3_atom y, u3_atom x) { _math_rnd=softfloat_round_near_even; return _rh_out(_rh_atan2(_rh_in(y), _rh_in(x))); }
  u3_noun u3wi_rh_atan2(u3_noun cor){ return _rh_jet2(cor, _rh_atan2); }
  u3_noun u3qi_rh_pow(u3_atom x, u3_atom n)   { _math_rnd=softfloat_round_near_even; return _rh_out(_rh_pow(_rh_in(x), _rh_in(n))); }
  u3_noun u3wi_rh_pow(u3_noun cor) { return _rh_jet2(cor, _rh_pow); }
  u3_noun u3qi_rh_pow_n(u3_atom x, u3_atom n) { _math_rnd=softfloat_round_near_even; return _rh_out(_rh_pow_n(_rh_in(x), _rh_in(n))); }
  u3_noun u3wi_rh_pow_n(u3_noun cor){ return _rh_jet2(cor, _rh_pow_n); }

#endif
