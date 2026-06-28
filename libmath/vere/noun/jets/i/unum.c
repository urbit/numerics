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

//  Elementary / transcendental functions.  Unary ones (incl. log-2/log-10,
//  which map to SoftUnum's base-2/base-10 log fns) reuse the unary template;
//  pow is binary; pow-n is posit ^ @u (the exponent is a raw integer, not a
//  posit).  Constants (pi/e/...) are left to pure Hoon, like /lib/math.
_UNUM_UNOP(exp, p8_exp, p16_exp, p32_exp)
_UNUM_UNOP(sin, p8_sin, p16_sin, p32_sin)
_UNUM_UNOP(cos, p8_cos, p16_cos, p32_cos)
_UNUM_UNOP(tan, p8_tan, p16_tan, p32_tan)
_UNUM_UNOP(log, p8_log, p16_log, p32_log)
_UNUM_UNOP(log2, p8_log2, p16_log2, p32_log2)
_UNUM_UNOP(log10, p8_log10, p16_log10, p32_log10)
_UNUM_UNOP(cbrt, p8_cbrt, p16_cbrt, p32_cbrt)
_UNUM_UNOP(atan, p8_atan, p16_atan, p32_atan)
_UNUM_UNOP(asin, p8_asin, p16_asin, p32_asin)
_UNUM_UNOP(acos, p8_acos, p16_acos, p32_acos)
_UNUM_UNOP(factorial, p8_factorial, p16_factorial, p32_factorial)
_UNUM_BINOP(pow, p8_pow, p16_pow, p32_pow)

/* ++pow-n:pp -- posit ^ @u (integer power); the exponent is a raw unsigned
** integer, not a posit, so it is read as a plain chub.
*/
  u3_noun
  u3qi_unum_pow_n(c3_d bloq, u3_atom x, u3_atom p)
  {
    c3_d ux = u3r_chub(0, x), up = u3r_chub(0, p), r;
    switch ( bloq ) {
      case 3:  r = p8_pow_n((posit8_t)ux, up);  break;
      case 4:  r = p16_pow_n((posit16_t)ux, up); break;
      case 5:  r = p32_pow_n((posit32_t)ux, up); break;
      default: return u3_none;
    }
    return u3i_chubs(1, &r);
  }

  u3_noun
  u3wi_unum_pow_n(u3_noun cor)
  {
    u3_noun x, p;  c3_d bloq;
    if ( c3n == u3r_mean(cor, u3x_sam_2, &x, u3x_sam_3, &p, 0) ||
         c3n == u3ud(x) || c3n == u3ud(p) ) return u3m_bail(c3__exit);
    if ( c3n == _unum_bloq(cor, &bloq) ) return u3_none;
    return u3qi_unum_pow_n(bloq, x, p);
  }

//  Rounding to integral value.  The Hoon rnd/flr/cel are eta-expanded into
//  unary gates wrapping +round, so they jet via the unary template.
_UNUM_UNOP(rnd, p8_nearest_int, p16_nearest_int, p32_nearest_int)
_UNUM_UNOP(flr, p8_floor, p16_floor, p32_floor)
_UNUM_UNOP(cel, p8_ceil, p16_ceil, p32_ceil)

/* ++sun:pp -- @u -> posit.  The argument is a raw unsigned integer (read as a
** full chub), NOT a posit pattern, so it is not masked to the posit width.
*/
  u3_noun
  u3qi_unum_sun(c3_d bloq, u3_atom v)
  {
    c3_d uv = u3r_chub(0, v), r;
    switch ( bloq ) {
      case 3:  r = p8_from_u64(uv);  break;
      case 4:  r = p16_from_u64(uv); break;
      case 5:  r = p32_from_u64(uv); break;
      default: return u3_none;
    }
    return u3i_chubs(1, &r);
  }
  u3_noun
  u3wi_unum_sun(u3_noun cor)
  {
    u3_noun v = u3r_at(u3x_sam, cor);  c3_d bloq;
    if ( u3_none == v || c3n == u3ud(v) ) return u3m_bail(c3__exit);
    if ( c3n == _unum_bloq(cor, &bloq) ) return u3_none;
    return u3qi_unum_sun(bloq, v);
  }

/* ++san:pp -- @s -> posit.  Decode the Hoon signed atom (even 2m -> +m,
** odd 2m-1 -> -m) to a C int64, then encode.
*/
  u3_noun
  u3qi_unum_san(c3_d bloq, u3_atom v)
  {
    c3_d  uv = u3r_chub(0, v);
    c3_ds sv = (uv & 1) ? -(c3_ds)((uv + 1) >> 1) : (c3_ds)(uv >> 1);
    c3_d  r;
    switch ( bloq ) {
      case 3:  r = p8_from_i64(sv);  break;
      case 4:  r = p16_from_i64(sv); break;
      case 5:  r = p32_from_i64(sv); break;
      default: return u3_none;
    }
    return u3i_chubs(1, &r);
  }
  u3_noun
  u3wi_unum_san(u3_noun cor)
  {
    u3_noun v = u3r_at(u3x_sam, cor);  c3_d bloq;
    if ( u3_none == v || c3n == u3ud(v) ) return u3m_bail(c3__exit);
    if ( c3n == _unum_bloq(cor, &bloq) ) return u3_none;
    return u3qi_unum_san(bloq, v);
  }

/* ++toi:pp -- posit -> (unit @s).  NaR -> ~ (none); else [~ @s] with the
** integer re-encoded into a Hoon signed atom (+m -> 2m, -m -> 2m-1).
*/
  u3_noun
  u3qi_unum_toi(c3_d bloq, u3_atom p)
  {
    c3_d  up = u3r_chub(0, p);
    c3_ds out;  c3_t ok;
    switch ( bloq ) {
      case 3:  ok = p8_to_i64((posit8_t)up, (int64_t*)&out);  break;
      case 4:  ok = p16_to_i64((posit16_t)up, (int64_t*)&out); break;
      case 5:  ok = p32_to_i64((posit32_t)up, (int64_t*)&out); break;
      default: return u3_none;
    }
    if ( !ok ) return u3_nul;
    c3_d sa = (out >= 0) ? ((c3_d)out << 1) : (((c3_d)(-out) << 1) - 1);
    return u3nc(u3_nul, u3i_chubs(1, &sa));
  }
  u3_noun
  u3wi_unum_toi(u3_noun cor)
  {
    u3_noun p = u3r_at(u3x_sam, cor);  c3_d bloq;
    if ( u3_none == p || c3n == u3ud(p) ) return u3m_bail(c3__exit);
    if ( c3n == _unum_bloq(cor, &bloq) ) return u3_none;
    return u3qi_unum_toi(bloq, p);
  }

/* ++is-close:pp -- |a - b| <= tol, a loobean.  Ternary sample [a b tol].
*/
  u3_noun
  u3qi_unum_is_close(c3_d bloq, u3_atom a, u3_atom b, u3_atom tol)
  {
    c3_d ua = u3r_chub(0, a), ub = u3r_chub(0, b), ut = u3r_chub(0, tol);  c3_t v;
    switch ( bloq ) {
      case 3:  v = p8_is_close((posit8_t)ua, (posit8_t)ub, (posit8_t)ut);    break;
      case 4:  v = p16_is_close((posit16_t)ua, (posit16_t)ub, (posit16_t)ut); break;
      case 5:  v = p32_is_close((posit32_t)ua, (posit32_t)ub, (posit32_t)ut); break;
      default: return u3_none;
    }
    return v ? c3y : c3n;
  }
  u3_noun
  u3wi_unum_is_close(u3_noun cor)
  {
    u3_noun a, b, tol;  c3_d bloq;
    if ( c3n == u3r_mean(cor, u3x_sam_2, &a, u3x_sam_6, &b, u3x_sam_7, &tol, 0) ||
         c3n == u3ud(a) || c3n == u3ud(b) || c3n == u3ud(tol) ) {
      return u3m_bail(c3__exit);
    }
    if ( c3n == _unum_bloq(cor, &bloq) ) return u3_none;
    return u3qi_unum_is_close(bloq, a, b, tol);
  }
