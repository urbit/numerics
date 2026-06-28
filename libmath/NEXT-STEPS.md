# `/lib/unum` — next steps

Status as of 2026-06-01.  The posit scalar layer (PR #13, merged) plus the
follow-up domain fixes (PR #38) implement a complete, SoftPosit-verified
`posit<n,2>` library: encode/decode, arithmetic, sqrt/round/fma, integer and
IEEE-754 conversion, the quire + fused dot product, and naive transcendentals,
at widths posit8/16/32/64/128.  This file tracks what is deliberately *not* yet
done, roughly in dependency order.

## Verification & tooling

- **Oracle**: SoftPosit (vendored C in `src/SoftPosit`, plus the `softposit`
  pip package) is the reference for arithmetic and quire ops.  It has **no
  transcendentals** — use `mpmath` → `convertDoubleToPX2` for those.
  `libmath/tools/posit_check.py` is the offline harness (exhaustive posit8,
  sampled posit16/32).
- **On-ship tests** run on a live ship (we used `~hex` via MCP).  There is no
  local urbit binary; a persistent build failure from the MCP tool is a *real*
  compile error — read the dojo (`tmux capture-pane`) for the actual
  `-need`/`-have`, since the MCP surface only returns a generic failure.
- A full per-width Hoon sweep of all 65k posit8 *pairs* times out the test
  runner; keep heavy exhaustion in the Python harness and keep the on-ship
  suite to round-trips, property checks, and curated spot values.

## 1. Standard-name alias layer  (small, additive)

The README's "Posit Standard Compliance" section is the spec: thin arms over
the implemented core.  Most are renames or one-line compositions:

- `negate`/`addition`/`subtraction`/`multiplication`/`division`,
  `compare-{equal,not-equal,greater,greater-equal,less,less-equal}`,
  `sign` ← `sgn`, `nearest-int`/`ceil`/`floor` ← `rnd`/`cel`/`flr`.
- `next`/`prior` — lexicographic successor/predecessor of the bit pattern
  (`+(p)` / `(dec p)` with NaR/extreme handling); genuinely new, but trivial.
- `*-pi` trig (`sin-pi` = `(sin (mul pi x))`), `*-plus-1`/`*-minus-1`
  elementary, `compound`/`root-n`, `hypot`, `fmm`, `arctan2` — compositions.
- inverse and hyperbolic trig (`arcsin`/`arccos`/`arctan`, `sinh`/`cosh`/`tanh`,
  `arcsinh`/…): new naive series, same caveats as the existing transcendentals.

No new infrastructure; can land as one PR.  Decide whether aliases live in
`unum.hoon` itself or a thin `unum-std.hoon` wrapper.

## 2. Transcendental accuracy  (medium)

The current `exp`/`sin`/`cos`/`log` are naive fixed-term Taylor series, accurate
only near the expansion point (documented).  To make them usable across the full
dynamic range:

- **Range reduction**: reduce `exp`/`log` by powers of 2 (posits make this
  exact), trig by multiples of `pi`/2.  Without it, `log(maxpos)` is off by ~3×.
- **Quire-accumulated sums**: run each series sum through the quire (`q-mul-add`
  → one `q-to-p`) so only the final rounding loses precision.  The hooks exist.
- Oracle: `mpmath` at high precision → round to the target posit.  Decide a
  per-width accuracy target (the 2022 standard requires correct rounding for
  compliance; we may settle for ≤1 ulp as an interim).

## 3. Lagoon `%unum` integration  (medium; the high-value item)

This is what makes posits useful for arrays/linear algebra.  Lagoon lives on the
**base desk** (`/lib/lagoon`, `/sur/lagoon`) — there is no separate lagoon desk.

- `sur/lagoon.hoon`: the `+$kind` union already has a commented-out `%unum`
  line — uncomment it.  Reconcile the stale aura note there (`@ruw/@ruh/@rub`)
  with the shipped family `@rpb/@rph/@rps/@rpd` (`bloq` selects width).
- `lib/lagoon.hoon` `fun-scalar`: add a `%unum` branch dispatching per-scalar to
  `add:rpb:unum`, `lth:rpb:unum`, etc., keyed by `bloq` (mirror the existing
  `%i754`/`%uint` branches).  Direct reuse of the implemented arms.
- **Reductions/linalg are the payoff**: route `dot`/`mmul`/`sum`/`cumsum`/
  `trace` over `%unum` arrays through the **quire** (`fdp` / repeated
  `q-mul-add` then one `q-to-p`) so products accumulate exactly and round once —
  exact dot product → matmul with no error accumulation.
- `convert`: posit↔i754 via the any-width `to-r*`/`from-r*` matrix (cross-width
  is fine and intended); posit↔uint/int2 via `sun`/`san`/`toi`.
- Jets come later (see §5).

## 4. Valids  (large; the third unum class)

The interval class (`@rvb/@rvh/@rvs`).  Entirely new surface: an interval is a
pair of posit endpoints with open/closed tags; arithmetic is interval
arithmetic.  Lowest priority; design from Gustafson's Type-III definition.

## 5. Jets  (separate, C/vere effort)

The library is pure Hoon by design.  SoftPosit's C (`src/SoftPosit/source/*.c`:
`pX2_*`, `qX2_*`, conversions) is the reference implementation and the spec for a
future jet, vendored into vere like SoftBLAS.  Only worth it once a consumer
(e.g. lagoon `%unum` matmul) needs the speed.

## Known minor items

- `fdp` silently truncates to the shorter of its two input lists.  Acceptable
  (zip semantics) but undocumented; add a note or `?>` if a caller needs strict
  equal-length.
- Cleanup the review suggested but we deferred: the eight `to-r*`/`from-r*`
  one-liners could collapse to a single `bloq`-dispatched gate; the per-width
  `test-consts-*` / round-trip tests could be table-driven (cf. the lagoon
  test-by-category convention).  Cosmetic.
- `twoc.hoon` now has corrected `lth`/`lte`/`gte` and `overflow`; the only
  former consumer (`lagoon-old.hoon`) was deleted in PR #38.  A real
  `%int2`-into-lagoon effort would be twoc's first live caller.

---

# Appendix: Valids (Type-III interval unums) — assessment & open questions

Status as of 2026-06-28.  Valids are the third Type-III unum (posits · quires ·
**valids**), the rigorous-interval class.  `/lib/unum` does not implement them;
this appendix records the design assessment so the eventual effort starts from a
shared understanding rather than a blank page.  Lowest priority (no consumer, no
oracle, loosely standardized) — slot it *after* the posit/quire jets.

## What a valid is

A valid is a *guaranteed enclosure* of a real quantity: two posit-like
**endpoints, each tagged with a "ubit" (uncertainty bit)**.  ubit = 0 → the
endpoint is exact (closed); ubit = 1 → "the open interval between this posit and
its neighbour".  So a valid encodes `[lo, hi]` / `(lo, hi)` / half-open — a value
*and* its uncertainty.  Every operation returns the **tightest valid containing
all possible results** (Gustafson's "end of error": a provable bound, not a
single rounded answer).

**The projective twist.**  Valids live on the *projective* real line — the reals
closed into a ring with a single point at the top (the NaR/∞ point, shared with
the posit bit layout).  This lets a valid represent **exterior / wrapping
intervals** ("x < −3 OR x > 5") and unbounded sets ("x > 5") by going the long
way round through infinity.  Consequence: compare / union / intersection are
arithmetic on *arcs of a circle*, not segments of a line — genuinely different
semantics from classical `[lo,hi]` interval arithmetic.

## What it would build on (and the one real new primitive)

Endpoints are posits, so valid arithmetic reuses the SoftUnum / `/lib/unum`
posit core.  The new building block is **directed rounding**: our posit `+bit`
only rounds nearest-even, but valid endpoints must round *outward* — the lower
bound toward −∞, the upper toward +∞ — to keep the enclosure sound.  So we need:

  - a round-toward-±∞ posit encoder (round-down / round-up of an exact value),
  - `+next` / `+prior` (lexicographic successor / predecessor of a posit bit
    pattern — already a TODO in §1 above; `+(p)` / `(dec p)` with NaR/extreme
    handling),
  - the ubit then records residual openness after directed rounding.

Everything else (endpoint add/sub/mul/div) is the existing posit arithmetic.

## Open design questions (pin these BEFORE writing code)

1. **Layout / width.**  Does `@rvb`/`@rvh`/`@rvs` mean a valid *built from two
   `posit{8,16,32}` endpoints*, or Gustafson's "valid⟨2n⟩ = two n-bit
   ubit-posits" packed into a byte/half/single?  The byte/half/single aura names
   mirror posits and are currently ambiguous.  **First decision.**
2. **Projective vs bounded.**  Full Gustafson (wrapping / exterior intervals on
   the projective ring) — the "real" valid — vs a simpler bounded `[lo,hi]`
   interval arithmetic that covers most numeric use with far less machinery.
3. **Set operations & special values.**  intersection / union / complement /
   `is-empty` / `is-everything`, plus the special encodings (empty set,
   all-reals, NaR).
4. **Comparisons.**  Posit ordering == two's-complement of the raw bits; that
   identity does **not** carry to valids — comparison becomes set / containment
   relations, needs its own design.

## Verification posture (harder than posits)

- **No SoftPosit oracle** — SoftPosit is posits + quires only.  Plan mirrors the
  transcendentals: a from-scratch **exact-rational interval reference** (Python
  intervals of `Fraction`s with directed rounding), and optionally **Stillwater
  `Universal`** (C++), which does implement `valid<>` types, as a second oracle.
- **Standardization caveat.**  The 2022 Posit Standard normatively pins posits
  and quires; valids are described in Gustafson's broader work but are *not*
  specified to the same degree.  So this is a clean-room *design* of both the
  Hoon `/lib/unum` arms and the C — NOT a transliteration of a fixed spec, which
  is a different (looser) posture from everything done so far, where `/lib/unum`
  was the bit-exact target.
