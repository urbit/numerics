// rd_exp_check.c -- SoftFloat-f64 reference for @rd exp (the jet body).
// Transcribed line-for-line from libmath/desk/lib/math.hoon ++rd ++exp
// (Cody-Waite reduction + degree-11 minimax Horner + scale2 ldexp).
// SoftFloat f64 round-near-even == the Hoon stdlib @rd ops, so this should be
// BIT-EXACT to (exp:rd:math x).  Prints "<in_hex> <out_hex>" per input.
//
// Build: re-archive the zig softfloat .a (not 8-byte aligned for Apple ld):
//   ar x <zig libsoftfloat.a> && libtool -static -o /tmp/libsoftfloat.a *.o
//   cc -I<sf-include> rd_exp_check.c /tmp/libsoftfloat.a -o /tmp/rd_exp_check
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include "softfloat.h"

typedef float64_t d;
static inline d   D(uint64_t b){ d r; r.v = b; return r; }
static inline d   mul(d a, d b){ return f64_mul(a,b); }
static inline d   sub(d a, d b){ return f64_sub(a,b); }
static inline d   add(d a, d b){ return f64_add(a,b); }

// pow2(j) = (j+1023)<<52 as f64 bits = 2^j  (normal range j in [-1022,1023])
static inline d pow2(int64_t j){ return D(((uint64_t)(j + 1023)) << 52); }

// scale2: correctly-rounded ldexp with overflow/subnormal tails (math.hoon:1519)
static d scale2(d p, int64_t k){
  if ( k - 1024 >= 0 )                       // k>=1024
    return mul(mul(p, pow2(1023)), pow2(k - 1023));
  if ( !(k + 1022 >= 0) )                     // k<-1022  (i.e. k<=-1023)
    return mul(mul(p, pow2(k + 54)), pow2(-54));
  return mul(p, pow2(k));
}

static d expd(d x){
  const uint64_t QNAN = 0x7ff8000000000000ULL, PINF = 0x7ff0000000000000ULL,
                 NINF = 0xfff0000000000000ULL;
  if ( !f64_eq(x, x) )       return D(QNAN);            // NaN
  if ( x.v == PINF )         return D(PINF);            // +inf
  if ( x.v == NINF )         return D(0);               // -inf -> 0
  d log2e = D(0x3ff71547652b82feULL);
  d ln2hi = D(0x3fe62e42fee00000ULL);
  d ln2lo = D(0x3dea39ef35793c76ULL);
  int64_t k = f64_to_i64(mul(x, log2e), softfloat_round_near_even, false);
  if ( k - 1025 >= 0 )       return D(PINF);            // overflow -> inf
  if ( !(k + 1075 >= 0) )    return D(0);               // underflow -> 0
  d ka = ui64_to_f64( (uint64_t)(k < 0 ? -k : k) );     // |k| as f64 (sun abs)
  d kf = (k >= 0) ? ka : sub(D(0), ka);
  d r  = sub( sub(x, mul(kf, ln2hi)), mul(kf, ln2lo) );
  // degree-11 minimax coeffs c0..c11 (math.hoon:1542-1547)
  static const uint64_t cs[12] = {
    0x3ff0000000000000ULL, 0x3ff0000000000000ULL, 0x3fe0000000000011ULL,
    0x3fc555555555555aULL, 0x3fa555555554f0cfULL, 0x3f8111111110f225ULL,
    0x3f56c16c187fbe02ULL, 0x3f2a01a01b14378fULL, 0x3efa01991ac8730aULL,
    0x3ec71ddf5749d126ULL, 0x3e928b4057f44145ULL, 0x3e5af631d0059becULL
  };
  // Horner: acc=0; for c in flop(cs) [c11..c0]: acc = acc*r + c   (math.hoon:1549)
  d p = D(0);
  for ( int i = 11; i >= 0; i-- ) p = add(mul(p, r), D(cs[i]));
  return scale2(p, k);
}

int main(int argc, char** argv){
  softfloat_roundingMode = softfloat_round_near_even;
  // inputs as native doubles (native == IEEE f64 bit patterns here)
  double xs[] = { 0.0, 0.5, 1.0, 2.0, 3.0, -1.0, -2.0, 5.0, 10.0, -10.0,
                  0.1, -0.1, 3.141592653589793, 100.0, -100.0,
                  700.0, -700.0, 709.0, -740.0 };
  int n = (int)(sizeof xs / sizeof xs[0]);
  for ( int i = 0; i < n; i++ ){
    union { double f; uint64_t b; } u; u.f = xs[i];
    d out = expd(D(u.b));
    printf("0x%016llx 0x%016llx  (x=%g)\n",
           (unsigned long long)u.b, (unsigned long long)out.v, xs[i]);
  }
  return 0;
}
