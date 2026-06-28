/// @file
///
/// Jets for the numerics `/lib/unum` library (2022 Posit Standard; userspace,
/// registered under the `non` chapter alongside `math` and `lagoon`).  Each jet
/// calls SoftUnum (ext/softunum), the bit-exact C twin of `/lib/unum`, so jet
/// output is identical to the unjetted Hoon.
///
/// `/lib/unum` is one generic posit core `++pp` (`|_ =bloq`) specialized into
/// `rpb`/`rph`/`rps`/`rpd`/`rpq` (posit8/16/32/64/128) via `%*`.  We register a
/// SINGLE `%unum` core and dispatch on `bloq` read from the door sample (the
/// same one-jet-per-op, runtime-dispatch shape `lagoon` uses).  bloq lives in
/// the `pp` door, which is the gate's context: gate axis 7 = door, door axis 6
/// = `=bloq`, so bloq is at gate axis 30.  SoftUnum covers bloq 3/4/5
/// (posit8/16/32); for bloq 6/7 (posit64/128, not yet in SoftUnum) the jet
/// returns u3_none and the pure-Hoon arm runs.
///
/// Marshalling uses chub reads/writes (word-size-agnostic across the 32- and
/// 64-bit runtimes).  Posit bit patterns occupy the low n bits.
///
/// MASTER COPY lives in urbit/numerics libmath/vere/noun/jets/i/unum.c; applied
/// by hand to the vere runtime.  SoftUnum itself is vendored (ext/softunum).

#include "jets/q.h"
#include "jets/w.h"
#include "noun.h"
#include "softunum.h"

//  bloq from the pp door sample: gate axis 7 = door, door axis 6 = =bloq.
#define _UNUM_BLOQ_AXIS 30

//  Read bloq (the door sample) from a gate core; c3n if absent/not-atom.
static inline c3_t
_unum_bloq(u3_noun cor, c3_d* out)
{
  u3_noun b = u3r_at(_UNUM_BLOQ_AXIS, cor);
  if ( u3_none == b || c3n == u3ud(b) ) {
    return c3n;
  }
  *out = u3r_chub(0, b);
  return c3y;
}

//  binary posit -> posit op (add/sub/mul/div), bloq-dispatched.
#define _UNUM_BINOP(nam, f8, f16, f32)                                       \
  u3_noun u3qi_unum_##nam(c3_d bloq, u3_atom a, u3_atom b) {                  \
    c3_d ua = u3r_chub(0, a), ub = u3r_chub(0, b), r;                         \
    switch ( bloq ) {                                                        \
      case 3:  r = f8((posit8_t)ua, (posit8_t)ub);    break;                 \
      case 4:  r = f16((posit16_t)ua, (posit16_t)ub); break;                 \
      case 5:  r = f32((posit32_t)ua, (posit32_t)ub); break;                 \
      default: return u3_none;                                               \
    }                                                                       \
    return u3i_chubs(1, &r);                                                 \
  }                                                                          \
  u3_noun u3wi_unum_##nam(u3_noun cor) {                                     \
    u3_noun a, b;  c3_d bloq;                                                \
    if ( c3n == u3r_mean(cor, u3x_sam_2, &a, u3x_sam_3, &b, 0) ||            \
         c3n == u3ud(a) || c3n == u3ud(b) ) return u3m_bail(c3__exit);       \
    if ( c3n == _unum_bloq(cor, &bloq) ) return u3_none;                     \
    return u3qi_unum_##nam(bloq, a, b);                                      \
  }

//  binary posit -> loobean comparison; returns & (c3y) / | (c3n).
#define _UNUM_CMP(nam, f8, f16, f32)                                         \
  u3_noun u3qi_unum_##nam(c3_d bloq, u3_atom a, u3_atom b) {                  \
    c3_d ua = u3r_chub(0, a), ub = u3r_chub(0, b);  c3_t v;                  \
    switch ( bloq ) {                                                        \
      case 3:  v = f8((posit8_t)ua, (posit8_t)ub);    break;                 \
      case 4:  v = f16((posit16_t)ua, (posit16_t)ub); break;                 \
      case 5:  v = f32((posit32_t)ua, (posit32_t)ub); break;                 \
      default: return u3_none;                                               \
    }                                                                       \
    return v ? c3y : c3n;                                                    \
  }                                                                          \
  u3_noun u3wi_unum_##nam(u3_noun cor) {                                     \
    u3_noun a, b;  c3_d bloq;                                                \
    if ( c3n == u3r_mean(cor, u3x_sam_2, &a, u3x_sam_3, &b, 0) ||            \
         c3n == u3ud(a) || c3n == u3ud(b) ) return u3m_bail(c3__exit);       \
    if ( c3n == _unum_bloq(cor, &bloq) ) return u3_none;                     \
    return u3qi_unum_##nam(bloq, a, b);                                      \
  }

//  unary posit -> posit op (neg/abs/sgn/sqt), bloq-dispatched.
#define _UNUM_UNOP(nam, f8, f16, f32)                                        \
  u3_noun u3qi_unum_##nam(c3_d bloq, u3_atom a) {                            \
    c3_d ua = u3r_chub(0, a), r;                                             \
    switch ( bloq ) {                                                        \
      case 3:  r = f8((posit8_t)ua);  break;                                 \
      case 4:  r = f16((posit16_t)ua); break;                                \
      case 5:  r = f32((posit32_t)ua); break;                                \
      default: return u3_none;                                               \
    }                                                                       \
    return u3i_chubs(1, &r);                                                 \
  }                                                                          \
  u3_noun u3wi_unum_##nam(u3_noun cor) {                                     \
    u3_noun a = u3r_at(u3x_sam, cor);  c3_d bloq;                            \
    if ( u3_none == a || c3n == u3ud(a) ) return u3m_bail(c3__exit);         \
    if ( c3n == _unum_bloq(cor, &bloq) ) return u3_none;                     \
    return u3qi_unum_##nam(bloq, a);                                         \
  }

_UNUM_BINOP(add, p8_add, p16_add, p32_add)
_UNUM_BINOP(sub, p8_sub, p16_sub, p32_sub)
_UNUM_BINOP(mul, p8_mul, p16_mul, p32_mul)
_UNUM_BINOP(div, p8_div, p16_div, p32_div)

_UNUM_CMP(lth, p8_lt, p16_lt, p32_lt)
_UNUM_CMP(lte, p8_le, p16_le, p32_le)
_UNUM_CMP(gth, p8_gt, p16_gt, p32_gt)
_UNUM_CMP(gte, p8_ge, p16_ge, p32_ge)
_UNUM_CMP(equ, p8_eq, p16_eq, p32_eq)
_UNUM_CMP(neq, !p8_eq, !p16_eq, !p32_eq)

_UNUM_UNOP(neg, p8_neg, p16_neg, p32_neg)
_UNUM_UNOP(abs, p8_abs, p16_abs, p32_abs)
_UNUM_UNOP(sgn, p8_sgn, p16_sgn, p32_sgn)
_UNUM_UNOP(sqt, p8_sqrt, p16_sqrt, p32_sqrt)

/* ++fma:pp -- fused multiply-add (a*b + c), single rounding.  Ternary gate:
** sample [a b c] = [a [b c]]: a @ sam_2, b @ sam_6, c @ sam_7.
*/
  u3_noun
  u3qi_unum_fma(c3_d bloq, u3_atom a, u3_atom b, u3_atom c)
  {
    c3_d ua = u3r_chub(0, a), ub = u3r_chub(0, b), uc = u3r_chub(0, c), r;
    switch ( bloq ) {
      case 3:  r = p8_fma((posit8_t)ua, (posit8_t)ub, (posit8_t)uc);    break;
      case 4:  r = p16_fma((posit16_t)ua, (posit16_t)ub, (posit16_t)uc); break;
      case 5:  r = p32_fma((posit32_t)ua, (posit32_t)ub, (posit32_t)uc); break;
      default: return u3_none;
    }
    return u3i_chubs(1, &r);
  }

  u3_noun
  u3wi_unum_fma(u3_noun cor)
  {
    u3_noun a, b, c;  c3_d bloq;
    if ( c3n == u3r_mean(cor, u3x_sam_2, &a, u3x_sam_6, &b, u3x_sam_7, &c, 0) ||
         c3n == u3ud(a) || c3n == u3ud(b) || c3n == u3ud(c) ) {
      return u3m_bail(c3__exit);
    }
    if ( c3n == _unum_bloq(cor, &bloq) ) return u3_none;
    return u3qi_unum_fma(bloq, a, b, c);
  }
