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

/* @rd log -- math.hoon ++rd ++log
**   x = 2^e * m, m in [sqrt(1/2),sqrt(2)); log(x) = e*ln2 + log(1+f),
**   f = m-1, s = f/(2+f), log(1+f) = f - s*(f - 2z*P2(z)), z=s*s, P2 the
**   atanh series 1/3 + z/5 + z^2/7 + ...
*/
  static float64_t
  _rd_log(float64_t x)
  {
    static const c3_d cs[10] = {       // atanh-series coeffs (math.hoon:1988)
      0x3fd5555555555555ULL, 0x3fc999999999999aULL, 0x3fc2492492492492ULL,
      0x3fbc71c71c71c71cULL, 0x3fb745d1745d1746ULL, 0x3fb3b13b13b13b14ULL,
      0x3fb1111111111111ULL, 0x3fae1e1e1e1e1e1eULL, 0x3faaf286bca1af28ULL,
      0x3fa8618618618618ULL
    };
    union doub r0, xx, m, f, s, z, p, c, rr, l1, efd, hi, lo, one, zero;

    zero.c = 0;
    r0.d = x;
    if ( !f64_eq(x, x) )                       { r0.c = _RD_QNAN; return r0.d; }  // NaN
    if ( r0.c == _RD_PINF )                    { return x; }                      // +inf
    if ( (r0.c == 0) || (r0.c == 0x8000000000000000ULL) ) { r0.c = _RD_NINF; return r0.d; }  // +-0 -> -inf
    if ( (r0.c >> 63) == 1 )                   { r0.c = _RD_QNAN; return r0.d; }  // x<0 -> NaN

    int    sub = ( ((r0.c >> 52) & 0x7ffULL) == 0 );
    xx = r0;
    if ( sub ) xx.d = f64_mul(x, _rd_bits(0x4350000000000000ULL).d);   // x * 2^54
    c3_ds  ae = sub ? -54 : 0;
    c3_d   b  = xx.c;
    c3_ds  ef = (c3_ds)((b >> 52) & 0x7ffULL) - 1023;
    m.c = (b & 0xfffffffffffffULL) | 0x3ff0000000000000ULL;            // m in [1,2)

    int big = f64_le(_rd_bits(0x3ff6a09e667f3bcdULL).d, m.d);          // m >= sqrt(2)
    if ( big ) { m.d = f64_mul(m.d, _rd_bits(0x3fe0000000000000ULL).d); ef = ef + 1; }
    ef = ef + ae;

    one.c = 0x3ff0000000000000ULL;
    f.d = f64_sub(m.d, one.d);
    s.d = f64_div(f.d, f64_add(m.d, one.d));
    z.d = f64_mul(s.d, s.d);

    p.c = 0;
    for ( int i = 10; i-- != 0; ) {            // Horner over flop(cs): c9..c0
      c.c = cs[i];
      p.d = f64_add(f64_mul(p.d, z.d), c.d);
    }
    rr.d = f64_mul(f64_add(z.d, z.d), p.d);                            // 2z * P2
    l1.d = f64_sub(f.d, f64_mul(s.d, f64_sub(f.d, rr.d)));             // f - s*(f - r)

    efd.d = ui64_to_f64( (c3_d)(ef < 0 ? -ef : ef) );
    if ( ef < 0 ) efd.d = f64_sub(zero.d, efd.d);                      // e as @rd
    hi.d = f64_mul(efd.d, _rd_bits(0x3fe62e42fee00000ULL).d);          // e*ln2hi
    lo.d = f64_mul(efd.d, _rd_bits(0x3dea39ef35793c76ULL).d);          // e*ln2lo
    return f64_add(hi.d, f64_add(l1.d, lo.d));
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

#endif
