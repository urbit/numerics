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

#include "jets/q.h"
#include "jets/w.h"

#include "noun.h"
#include "softfloat.h"

  union doub {
    float64_t d;
    c3_d c;
  };

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
    if ( !f64_eq(x, x) )                  { r0.c = 0x7ff8000000000000ULL; return r0.d; }  // NaN
    if ( r0.c == 0x7ff0000000000000ULL )  { return x; }                                   // +inf
    if ( r0.c == 0xfff0000000000000ULL )  { r0.c = 0; return r0.d; }                       // -inf -> 0

    log2e.c = 0x3ff71547652b82feULL;
    ln2hi.c = 0x3fe62e42fee00000ULL;
    ln2lo.c = 0x3dea39ef35793c76ULL;

    c3_ds k = (c3_ds)f64_to_i64(f64_mul(x, log2e.d), softfloat_round_near_even, 0);
    if ( (k - 1025) >= 0 )    { r0.c = 0x7ff0000000000000ULL; return r0.d; }  // overflow -> inf
    if ( !((k + 1075) >= 0) ) { r0.c = 0; return r0.d; }                       // underflow -> 0

    ka.d = ui64_to_f64( (c3_d)(k < 0 ? -k : k) );
    kf.d = (k >= 0) ? ka.d : f64_sub((union doub){.c=0}.d, ka.d);
    rr.d = f64_sub( f64_sub(x, f64_mul(kf.d, ln2hi.d)), f64_mul(kf.d, ln2lo.d) );

    p.c = 0;
    for ( c3_w i = 12; i-- != 0; ) {       // Horner over flop(cs): c11..c0
      c.c = cs[i];
      p.d = f64_add(f64_mul(p.d, rr.d), c.d);
    }
    return _rd_scale2(p.d, k);
  }

  u3_noun
  u3qi_rd_exp(u3_atom a)
  {
    union doub c, e;
    softfloat_roundingMode = softfloat_round_near_even;
    c.c = u3r_chub(0, a);
    e.d = _rd_exp(c.d);
    return u3i_chubs(1, &e.c);
  }

  u3_noun
  u3wi_rd_exp(u3_noun cor)
  {
    u3_noun x = u3r_at(u3x_sam, cor);

    if ( u3_none == x || c3n == u3ud(x) ) {
      return u3m_bail(c3__exit);
    }
    return u3qi_rd_exp(x);
  }
