# `/lib/rand` — next steps

Status as of 2026-07-12. Milestones 1-6 (`rand-spec.md` section 13) are done:
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
   - `/lib/complexrand` -- `+cuniform`, `+normal-parts`, `+cnormal`,
     `+on-circle`, `+in-disk`. Every arm here draws floats directly, so
     this one DOES depend on `/lib/i754rand` (for `+rd:uni`/`+rd-oo:uni`),
     plus `/lib/complex` and `/lib/math` (cos/sin).
   - `/lib/unumrand` -- `+posit-lattice`, `+posit-unit`. Needs only
     `/lib/rand` + `/lib/unum`; no floats (`+posit-unit`'s construction is
     explicitly float-free per rand-spec.md section 12.4).
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
     twocrand and fixedrand are DONE, complexrand remains); then unumrand
     alone, since `++posit-unit`'s value-uniformity claim needs a
     NEW mpmath oracle script (`librand/tools/posit_unit_check.py`,
     alongside `unum_cheb_check.py`) before it can be trusted, unlike
     everything shipped so far which either had a real external KAT or was
     ship-verified against a hand/Python-traced computation.
8. Saloon `+rand-ray` (counter-window design) + replay-equality tests.
9. Full README rewrite + this file's final pass (ziggurat, BTPE, buffered
   Philox, posit rays, quire Monte Carlo note, `@rh`/`@rq` distributions —
   all deliberately deferred out of v1, per `rand-spec.md`).
