# `/lib/rand` — next steps

Status as of 2026-07-13: **all nine milestones of `rand-spec.md` section 13
are done.** This file is a development log (decisions, footguns, bugs
found and fixed) kept in full rather than rewritten, per this project's own
convention of not silently erasing history — see "Deferred to v2" at the
bottom for what's actually left to do. Milestones 1-6 (`rand-spec.md`
section 13) landed first:
`++split-mix`, `++philox`, `++seed`, `+step`, `+fork`, `++gen`, `++uni`,
`++pcg`, `++dist` (both `++rd` and `++rs`), `++sample` -- KAT-verified
against Vigna's reference SplitMix64, the Random123 `kat_vectors` file, an
independent Python IEEE-754 encoder for `++uni`'s float arms, pcg-c's own
seed=42/seq=54 demo convention for `++pcg`, NumPy's `random_poisson_ptrs`
for `++dist`'s PTRS branch, and on-ship-computed regression values plus
moment tests elsewhere.

## Architecture: split into `/lib/rand` (plumbing) + `/lib/i754rand` (porcelain)

After milestone 6 landed, the single `rand.hoon` file was split in two,
matching how `/lib/twoc`/`fixed`/`complex`/`unum` are already kept separate
from `/lib/math` in this codebase:

- **`/sur/rand`** -- `+$phil`/`+$sm64`/`+$pcg64`/`+$rng`, extracted so any
  file needing just the type can use a lightweight `/-` import.
- **`/lib/rand`** (the plumbing) -- `++philox`, `++split-mix`, `++pcg`,
  `++seed`, `+step`, `+fork`, `++gen`, `++uni` (integer-only: `+bits`/
  `+below`/`+between`), `++sample` (shuffle/permutation/choice/reservoir).
  Every planned adapter library (`i754rand`, and the non-float ones below)
  imports this file. Tests: `-test %/tests/lib/rand ~` (44 tests).
- **`/lib/i754rand`** (the IEEE-754 porcelain) -- `++uni`'s float arms
  (`+rs`/`+rd`/`+rh`/`+rq`/`+rs-oo`/`+rd-oo`), `++alias` (Vose's method,
  moved here from `++sample`), and `++dist`. Tests:
  `-test %/tests/lib/i754rand ~` (47 tests).

`++alias` moved out of `++sample` specifically because `+draw` needs an
actual uniform `@rd` float draw to decide accept-vs-redirect -- it isn't
float-free the way the rest of `++sample` is, so it belongs with the
porcelain, not the plumbing. `+categorical` (in `++dist`) references it as
a same-file sibling (`alias-table`/`draw:alias`, unqualified).

Import mechanics worth remembering if this gets refactored again: `/+ rand`
alone does NOT make the unwrapped `/sur/rand` types resolve as `rng:rand`
externally unless `/lib/rand` itself does `=+ rand` after its own `/- rand`
import (confirmed empirically -- `rng:rand` works from any consuming file
once `rand.hoon` does the unwrap). Adding a second `=+ rand` inside a
*consumer* file (e.g. `i754rand.hoon`) to get `rng` bare, instead of
`rng:rand`, breaks qualified sibling lookups like `bits:uni:rand` in ways
that were not obvious from the error message (`-find.uni`) -- so
`i754rand.hoon` does NOT unwrap; every bare `rng` in its arm signatures is
written `rng:rand` explicitly instead.

**Spec correction found during milestone 4**: `rand-spec.md` section 3.3
originally said PCG64 outputs the current state's permutation and THEN
advances. That's backwards -- pcg-c's actual reference
(`pcg_setseq_128_xsl_rr_64_random_r`, cross-checked against NumPy's vendored
copy) advances the state FIRST and outputs from the new state. Both the
spec and `++pcg` now follow the verified reference order; see the struck-
through note left in the spec for the record.

**Two real algorithm bugs caught during milestone 6's on-ship testing**
(both are exactly why this library tests against ship-computed values,
not hand-derivation): `++alias`'s `+build-loop` had the small/large
worklist assignment reversed (an index whose weight said it belonged on
one worklist was pushed onto the other); `++dist`'s `+binomial` compared
the uniform draw against the individual pmf term instead of the
accumulated CDF, which either returned 0 too often or (for larger n) ran
the recurrence past n and crashed with subtract-underflow. Both are now
documented in the arms' own doccords so a future reader sees why the code
looks the way it does.

Note on `++dist`: `+chi2`/`+student-t` aren't in milestone 5's literal arm
list but are one-line compositions of `+gamma`/`+normal` with no new
machinery, so they were added alongside rather than left in limbo until an
unscheduled milestone. `+categorical` has no `++rs` mirror: alias-table's
`.prob` is fixed at `@rd`, so it's precision-invariant.

Remaining milestones, in dependency order (see `rand-spec.md` section 13 for
full detail):

7. Non-float adapters, as four SEPARATE libraries (not nested in `/lib/rand`
   or `/lib/i754rand`), each importing `/lib/rand` for the engine/uni/sample
   primitives:
   - DONE: `/lib/twocrand` -- `+twoc-full`, `+twoc-between`. Needs only
     `/lib/rand` + `/lib/twoc`; no floats at all.
   - DONE: `/lib/fixedrand` -- `+fixed`, `+fixed-unit`, `+fixed-between`
     (delegates to `twocrand`'s `+twoc-between`). Needs `/lib/rand` +
     `twocrand` + `/lib/fixed` (needed a `+from-rd` added, mirroring its
     existing `+from-rs` -- done). Footgun hit and documented: `+fixed`
     (the arm, matching rand-spec.md's public API) collides with `fixed`
     (the imported library face), since a same-named battery arm shadows
     an imported face throughout its own core -- same class of bug
     `/lib/fixed` itself hit renaming `+twoc` to `+neg`. Fixed here by
     aliasing the import (`/+ fx=fixed`) instead of renaming the arm,
     since the arm name is spec-mandated and the import name isn't.
   - DONE: `/lib/complexrand` -- `+cuniform`, `+normal-parts`, `+cnormal`,
     `+on-circle`, `+in-disk`. Every arm here draws floats directly, so
     this one DOES depend on `/lib/i754rand` (for `+rd:uni`/`+rd-oo:uni`),
     plus `/lib/complex` and `/lib/math` (cos/sin). One arm set per
     component-width door (`+cd` first, `+cs` mirror second, matching
     `/lib/complex`'s own ship order). Caught a real, previously-
     uncovered bug while building `+cnormal:cs`: `/lib/math`'s `@rs`
     `+invsqt2` was written `.70710677` (missing the leading `.0.`),
     which Hoon parses as the integer 70,710,677.0, not 0.70710677 --
     every other precision (`@rd`/`@rh`/`@rq`) had the correct form.
     Fixed, with a new regression suite (`tests/lib/math-constants.hoon`
     in `libmath`) covering `tau`/`pi`/`phi`/`sqt2`/`invsqt2` at all four
     precisions, since nothing previously exercised any of them.
   - DONE: `/lib/unumrand` -- `+posit-lattice`, `+posit-unit`. Needs only
     `/lib/rand` + `/lib/unum`; no floats (`+posit-unit`'s construction is
     explicitly float-free per rand-spec.md section 12.4). Distributions
     ("sample at @rd, convert") need NO new /lib/unum plumbing, unlike
     fixedrand's `+from-rd` -- `/lib/unum` already ships `+from-rh/rs/rd/rq`
     at every width door.
   - Range subtlety (worked through with the user before implementation):
     `+posit-lattice` (uniform over bit patterns, minus NaR) and
     `+posit-unit` (uniform over VALUES on [0,1), exact) are fundamentally
     different because posits are tapered -- consecutive bit patterns are
     NOT evenly spaced in value. `+posit-lattice` ships at all FIVE width
     doors (rpb/rph/rps/rpd/rpq); it has no bit-count subtlety (simple
     NaR-rejection, exact at any width).
   - `+posit-unit`'s bit-count k is the dangerous part: draw k raw bits u,
     encode the dyadic u*2^-k via /lib/unum's existing +bit (RNE,
     saturating -- confirmed no rounding-mode parameter needed). For this
     to be EXACTLY the round-to-nearest image of continuous uniform (not
     approximately), k must exceed the finest rounding-cell width anywhere
     in [0,1) -- which, because of the taper, is NEAR ZERO and shrinks fast
     with width. Hand-traced /lib/unum's own `+sea` decode of posit8's
     pattern `1` (minpos) and confirmed minpos = 2^-4(n-2) exactly (2^-24
     for n=8), so the rounding boundary nearest zero sits at 2^-(4(n-2)+1).
     The spec's k values (32/64/128 for posit8/16/32) are exactly k=4n,
     which gives a CONSTANT 7-bit safety margin at any width
     (4n - (4(n-2)+1) = 7 always) -- naive choices like k=n or k=n+8 would
     be silently wrong (they'd truncate resolution near zero and bias the
     distribution with no obvious symptom). u=0 must produce `[%z ~]`
     explicitly per spec (not a degenerate `[%p ...]` with a=0). Scoped to
     `[0,1)` only -- a general bounded-range posit uniform is out of scope.
   - Decided with the user: `+posit-unit` stays scoped to posit8/16/32
     (matching /lib/unum's existing "unverified until oracle sweep extends"
     caveat for posit64/128 -- even though the k=4n margin argument would
     work mathematically at any width, since +bit's encoding logic isn't a
     convergence-based transcendental, staying scoped avoids a false sense
     of rigor where the rest of the library hasn't independently verified
     those widths either).
   - Decided with the user on oracle rigor: chi-square posit8 (256 patterns,
     exactly enumerable) AND posit16 (65,536 patterns, still cheap to
     enumerate exactly) against the mpmath oracle; posit32 (2^32 patterns,
     exact enumeration infeasible) gets ship-verified value-regression
     tests only, relying on the same proven k=4n margin argument rather
     than independent re-verification.
   - Plan revised during implementation: originally planned as a real Hoon
     `-test` regression (hardcoded expected statistic/threshold). Dropped
     in favor of `librand/tools/posit_unit_check.py` staying a standalone
     oracle (no mpmath needed -- everything here is already exactly dyadic,
     so plain `Fraction` arithmetic is exact throughout): embedding a
     256- or 65536-row expected-probability table as Hoon literals isn't
     practical, and unlike the moment tests (which check a couple of
     summary statistics), this needs the FULL per-pattern distribution to
     mean anything. `posit_unit_check.py table <width>` prints the exact
     table; `posit_unit_check.py chi2 <width> <counts-file> --bins N`
     chi-squares ship-drawn raw patterns (one per line) against it,
     quantile-binned so N stays modest. The oracle script's OWN
     correctness is validated in-repo (a simulated correct-distribution
     sample passes with a sane p-value; a deliberately-wrong one is
     rejected with p~0 -- proving the test has real power, not just
     rubber-stamping). Actually run once at posit8 with N=100,000 ship-
     drawn draws (fixed seed 12345), tallied on-ship into a histogram (at
     most 65 possible patterns for posit8, so the histogram itself is
     small/printable) via a throwaway `.hoon` file (NOT a raw dojo one-
     liner -- a long single-line multi-`=/`/`|-` dojo command got silently
     truncated in transmission mid-session and left the dojo edit buffer
     stuck with an unclosed expression; recovered with Ctrl-E then Ctrl-U,
     not Ctrl-U alone, since the cursor was stuck at the line's start;
     lesson: bulk/long computations belong in a deployed `.hoon` file run
     via `-build-file`, not a giant single dojo line). Result: 61/65
     patterns hit (the 4 misses are the lowest-probability patterns near
     zero, expected counts well under 1 at this N), binned into 20
     quantile bins (some of the highest-probability patterns near 1.0
     individually exceed 1/32 of the mass, so `--bins 32` yields fewer
     than 32 actual bins -- expected, see the script's own comment),
     chi-square = 26.86 on 19 dof, p = 0.108 -- consistent with the exact
     expected distribution at any standard significance level. posit16's
     command is documented and ready to run the same way but wasn't
     executed in this pass (extracting a large enough on-ship sample is
     more awkward at that width, and its 65,536-entry histogram would be
     unwieldy to print/parse from dojo); the bit-exact cross-check against
     `posit_unit_check.py`'s `encode()` at posit16/32 (done for several
     seeds, in `+test-posit-unit`) is a stronger per-draw guarantee than
     the statistical test anyway, since it's the same deterministic code
     path at every width. To rerun or extend: deploy a `.hoon` file like
     `=/  count  N  =/  r  (from-atom:seed:rand %sm64 SEED)  =/  m
     *(map @ @)  =/  i  0  |-  ^-  (map @ @)  ?:  =(i count)  m  =^  v  r
     (posit-unit:rpb:unumrand r)  =/  c  (fall (~(get by m) v) 0)
     $(i +(i), m (~(put by m) v +(c)))`, `-build-file` it, `~(tap by
     <face>)` to print, paste the pairs into a Python dict, and feed to
     `expected_probabilities`/`quantile_bins`/`chi_square_stat` in
     `posit_unit_check.py` directly.
   - Decided: the "sample at @rd, quantize/convert" pattern mentioned in
     rand-spec.md 12.1/12.4 (fixed-point and posit *distributions*) is
     NOT shipped as dedicated per-distribution wrapper arms (that would be
     ~20 nearly-identical one-liners across fixedrand/unumrand doing
     nothing beyond composing two already-existing arms). Instead: document
     the one-line composition (`i754rand`'s dist arm -> the target type's
     own `+from-rd`) in each adapter's own doccord, and add ONE test per
     adapter proving the composition actually works (this also validates
     `+from-rd` once it's added to `/lib/fixed`). Revisit if real callers
     end up wanting the same composition in more than one or two places.
   - Split into two passes: twocrand + fixedrand + complexrand first
     (mechanical, same ship-verification rhythm as prior milestones --
     ALL THREE ARE DONE); then unumrand alone, since `++posit-unit`'s
     value-uniformity claim needs a
     NEW mpmath oracle script (`librand/tools/posit_unit_check.py`,
     alongside `unum_cheb_check.py`) before it can be trusted, unlike
     everything shipped so far which either had a real external KAT or was
     ship-verified against a hand/Python-traced computation.
8. DONE: Saloon `+rand-ray` (counter-window design) + replay-equality tests.
   Lands directly in `/lib/saloon`'s `+sa` core (a new `+|  %rand` section),
   not a separate librand file -- the spec's own framing ("In /lib/saloon,
   a core taking a Lagoon meta and an rng") describes per-arm shape, not a
   literal nested sub-core; Saloon's existing convention is one flat `+sa`
   door with `+|` section markers, so `+rand-ray` follows that rather than
   introducing complexrand/unumrand-style nesting.
   - `+fill-uniform` (single non-rejecting draw/element, %phil gets
     ctr0+i directly), `+fill-normal`/`+fill-expon` (rejection-based,
     %phil gets the ctr0+i*2^32 window with an explicit crash on
     exhaustion), `+fill-below` (%uint rays via Lemire, also windowed).
     `%i754` only, bloq 5/6 (`@rs`/`@rd`) -- matches rand-spec.md's stated
     v1 scope ("posit rays deferred").
   - Blocked on a pre-existing, unrelated bug: `+sa`'s scalar transcendental
     dispatch (`+trans-scalar`, `+fadd`/etc) called `/lib/math`'s doors as
     `[rnd rtol]` (2-tuple) where they need `[r rtol atol]` (3-tuple) --
     `saloon.hoon` didn't `-build-file` AT ALL before this, on ANY branch,
     confirmed by testing the unmodified file directly. Fixed and shipped
     as its OWN PR (#77, `sigilante/saloon-scalar-rtol-fix` off `main`),
     since it's also needed upstream in `urbit/urbit` independent of
     librand. That fix is ALSO carried on this branch (duplicated, not
     rebased) so `+rand-ray` has something to build against; reconcile via
     rebase once #77 merges to main.
   - Two Hoon footguns hit writing `+rand-ray` itself, both the "narrowing
     doesn't survive a recursive `$(...)` rebind" class already known in
     this codebase: (a) mutating `.r` (`r(ctr.p ...)`) inside a `?=(%phil
     -.r)`-narrowed recursive trap lost the narrowing on recursive re-
     entry -- fixed by building FRESH `[%phil key0 ctrN]` literals instead
     of mutating; (b) the non-%phil ("sequential engine") branch's `.r`
     is narrowed to EXCLUDE %phil by the same `?:`, but `+draw`'s return
     type is the full `rng` union, so the trap's first entry (narrow) and
     recursive re-entries (widened by `=^ v r (draw r)`) disagreed --
     fixed by explicitly widening `` `rng:rand`r `` once before the trap
     so every entry matches.
   - Every arm bit-exact cross-checked against a DIRECT call to the
     underlying `/lib/rand`/`/lib/i754rand` primitive at the expected
     counter (not just "doesn't crash") -- see `tests/lib/saloon-rand-
     ray.hoon` in the `saloon` desk (not `librand`, since the code lives
     in Saloon).
9. DONE: full README rewrite + this file's final pass (this section).

## Deferred to v2

Six items, all deliberately out of v1 scope per `rand-spec.md`. None of
these are bugs or gaps in what's shipped — they're follow-on work a future
session can pick up independently, in no particular order:

- **Ziggurat.** `+normal`/`+expon` use Marsaglia polar / inversion, which
  are simple and correctly rounded but do a `+sqt`/`+log` call per
  accepted sample (Marsaglia's rejection rate is ~21.5%, so ~1.27 draws
  and one transcendental call per deviate on average). Ziggurat avoids
  the transcendental call entirely in the common case via precomputed
  layer tables, at the cost of real implementation complexity (building
  and validating the tables) for a speedup that mostly matters once these
  arms are jetted. Not worth it before then.
- **BTPE for `+binomial`.** The current inversion-by-CDF-accumulation
  method crashes above `n*min(p,1-p) >= 30` (the recurrence gets slow and
  numerically risky past that point). Kachitvichyanukul & Schmeiser's
  BTPE algorithm (Binomial, Triangle, Parallelogram, Exponential regions)
  handles the large-n regime in O(1) expected time, but it's a
  substantially more involved algorithm than anything else in `++dist` —
  a real follow-on project, not a quick add.
- **Buffered Philox.** `+step`'s `%phil` branch keeps only the low 64
  bits of each 128-bit Philox4x32-10 block and discards the high 64 bits
  (`/lib/rand`'s own `+step` doc comment: "wasting them is acceptable at
  v1... a buffered variant that reuses both halves is a NEXT-STEPS
  item"). A buffered variant would cache the unused half and serve it on
  the *next* `+step` call instead of computing a fresh block, roughly
  doubling sequential throughput with zero extra Philox rounds. Doesn't
  change any output bit-for-bit; purely a performance follow-up, and one
  that interacts with jetting (the cache would need to live in the `rng`
  noun itself, changing `+$phil`'s shape) so it's more natural to do
  alongside milestone 11's jetting work than before it.
- **Posit rays.** Saloon's `+rand-ray` (milestone 8) only fills `%i754`
  rays. Extending `+fill-uniform`/`+fill-normal`/etc. to `%unum`-kind rays
  (posit-valued Lagoon arrays), using `/lib/unumrand`'s `+posit-unit`/
  `+posit-lattice` as the per-element generator, is mechanically similar
  to the existing `%i754` path but was out of scope for the milestone
  that shipped it.
- **Quire Monte Carlo.** Noted when `/lib/unumrand` shipped: `/lib/unum`'s
  quire (`+fdp`, Type III unums' exact fixed-point accumulator) makes
  sample *sums* singly-rounded — a capability hardware IEEE floats simply
  don't have (every float addition in an accumulation loop rounds; a
  quire-based reduction rounds exactly once, at the very end). This is a
  genuinely novel capability worth designing a Monte-Carlo-reduction API
  around, but it's a design project in its own right, not a small
  addition to an existing arm.
- **`@rh`/`@rq` distributions.** `/lib/i754rand`'s `++dist` only has
  `++rd` (reference, double)/`++rs` (single) mirrors, matching
  rand-spec.md's own scoping ("each arm exists at @rd (reference) and
  @rs"). Half (`@rh`) and quad (`@rq`) precision distributions were never
  in v1 scope; adding them is a mechanical re-instantiation of the same
  algorithms (matching how `++rs` itself was built as a mirror of
  `++rd`), gated only by whether a real caller needs sub-single or
  above-double precision sampling.
