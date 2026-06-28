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

/* ++add:pp -- posit add, dispatched on bloq (posit8/16/32 jetted).
*/
  u3_noun
  u3qi_unum_add(c3_d bloq, u3_atom a, u3_atom b)
  {
    c3_d ua = u3r_chub(0, a);
    c3_d ub = u3r_chub(0, b);
    c3_d r;
    switch ( bloq ) {
      case 3:  r = p8_add((posit8_t)ua, (posit8_t)ub);    break;
      case 4:  r = p16_add((posit16_t)ua, (posit16_t)ub); break;
      case 5:  r = p32_add((posit32_t)ua, (posit32_t)ub); break;
      default: return u3_none;                  //  posit64/128: fall back to Hoon
    }
    return u3i_chubs(1, &r);
  }

  u3_noun
  u3wi_unum_add(u3_noun cor)
  {
    u3_noun a, b, blq;
    if ( c3n == u3r_mean(cor, u3x_sam_2, &a, u3x_sam_3, &b, 0) ||
         c3n == u3ud(a) || c3n == u3ud(b) ) {
      return u3m_bail(c3__exit);
    }
    blq = u3r_at(_UNUM_BLOQ_AXIS, cor);
    if ( u3_none == blq || c3n == u3ud(blq) ) {
      return u3_none;
    }
    return u3qi_unum_add(u3r_chub(0, blq), a, b);
  }
