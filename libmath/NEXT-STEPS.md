# `/lib/unum` — next steps

Status as of 2026-06-01.  The posit scalar layer (PR #13, merged) plus the
follow-up domain fixes (PR #38) implement a complete, SoftPosit-verified
`posit<n,2>` library: encode/decode, arithmetic, sqrt/round/fma, integer and
IEEE-754 conversion, the quire + fused dot product, and naive transcendentals,
at widths posit8/16/32/64/128.  This file tracks what is deliberately *not* yet
done, roughly in dependency order.

## Roadmap — 2026-06-28 (post-jets)

The jet effort (old §5) is **done**: **SoftUnum** (`sigilante/SoftUnum`, the
pure-integer C twin of `/lib/unum`, verified four ways) is vendored into vere
(`ext/softunum`) and the full `/lib/unum` surface (54 arms, posit8/16/32) is
jetted under `non/unum`, fired bit-exact on a hoon-135 fakezod, and benchmarked
(arithmetic ~180–295×, transcendentals ~760–2000×, `fdp` ~760–2300× over
interpreted).  PRs: numerics #65 (hints + reference + valids appendix +
benchmarks); vere #1046 (draft, stacked on the math-jets PR #1044).

Prioritized next steps:

0.  **Land the PRs.**  numerics #65 → merge.  vere #1046 is gated on #1044
    (math jets) and #1022 (lagoon SoftBLAS) merging; then re-target #1046 to
    `develop` and mark ready.  Registration is hoon-135 only (lowest kelvin).

1.  **Posit linear algebra — the Saloon/Lagoon payoff (highest value).**
    `%unum` is already wired into Lagoon (PR #42) and Saloon eig (PR #58), so
    Saloon-over-posits works and now inherits the *scalar* jet speedup (its
    inner `mmul`/`dot`/`sqt`/`add` hit the jetted unum ops).  The remaining win
    is **array-level**: a SoftUnum-backed **posit GEMM/dot using the quire**
    (the analog of SoftBLAS), vendored like SoftBLAS, with a Lagoon `%unum`
    `dot`/`mmul` jet dispatching to it -- one exact-accumulated C call instead
    of a Hoon loop over per-element `fdp`.  This is the headline posit feature:
    exact dot product -> matmul with no error accumulation.  Saloon
    decompositions inherit it for free.

2.  **Perf gaps the benchmark surfaced.**
    - posit16/32 `sqt`/`atan` are slow even jetted (jetted `atan:rps` ~850 us):
      the 512-bit `wide_t` **bit-by-bit `isqt`** in the AGM loop is the
      bottleneck.  Replace with a faster wide sqrt (Newton + wide divmod, or an
      `__int128` seed then refine).
    - posit16 over-uses the 512-bit `wide_t` (its quire is 256-bit, arithmetic
      fits `__int128`); a tighter p16 path would shave the common arithmetic.

3.  **Coverage: posit64/128.**  SoftUnum does 8/16/32; `rpd`/`rpq` jets return
    `u3_none` (fall back to Hoon).  The `wide_t` machinery already exists --
    extend SoftUnum to 64/128 (1024/2048-bit quire) to jet the last two widths.
    No external oracle (cerlane `pX2` caps at 32); verify vs the from-scratch
    reference + the Hoon.

4.  **Accuracy: range-reduced / quire-accumulated transcendentals** (old §2).
    The naive Taylor `exp`/`log`/`sin` are accurate only near the expansion
    point.  Range-reduce (powers of 2 are exact for posits; trig by pi/2) and
    run the series through the quire.  Coordinated Hoon + SoftUnum change (both
    stay bit-exact).

5.  **Standard-name alias layer** (old §1): the 2022-standard public names
    (`addition`/`subtraction`/`sin-pi`/`compound`/`hypot`/`arctan2`/...) over
    the implemented core; mostly renames + thin compositions.

6.  **Valids** (the Type-III interval class) -- see the appendix below.  Lowest
    priority: new surface, no SoftPosit oracle, loosely standardized.

7.  **Parallel Rust SoftUnum** (for NockApp; not urgent).  C and Rust can't
    share code, so "same behavior" is enforced by a language-neutral
    test-vector corpus generated from the oracle that both reproduce.  The C
    repo reserves a top-level `rust/` seam.

8.  **Upstream SoftUnum -> `urbit/SoftUnum`** (like `urbit/SoftBLAS`) once it
    stabilizes, and re-pin the vere `ext/softunum` tarball there.

9.  **Decimal printer / literal syntax** for the posit auras (`@rpb`..`@rpq`):
    emit the §6.3 minimum significant digits.  Optional local runtime patch,
    decoupled from the library.

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

## 5. Jets  (DONE -- see the 2026-06-28 roadmap above)

**Done.**  Rather than vendor SoftPosit (whose `pX2` caps at 32 bits and has no
transcendentals), we built **SoftUnum** (`sigilante/SoftUnum`) -- the bit-exact
pure-integer C twin of this library -- vendored it into vere (`ext/softunum`),
and jetted the full surface (54 arms, posit8/16/32) under `non/unum`.  Verified
firing bit-exact on a hoon-135 fakezod; benchmarked in `benchmark/results/`.
PRs: numerics #65, vere #1046 (draft on #1044).  Follow-on perf/coverage work
(posit GEMM, wide `isqt`, posit64/128) is in the roadmap above.

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
