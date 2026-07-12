# `/lib/rand` — next steps

Status as of 2026-07-12. Milestones 1-6 (`rand-spec.md` section 13) are done:
`++split-mix`, `++philox`, `++seed`, `+step`, `+fork`, `++gen`, `++uni`,
`++pcg`, `++dist` (both `++rd` and `++rs`), `++sample` -- KAT-verified
against Vigna's reference SplitMix64, the Random123 `kat_vectors` file, an
independent Python IEEE-754 encoder for `++uni`'s float arms, pcg-c's own
seed=42/seq=54 demo convention for `++pcg`, NumPy's `random_poisson_ptrs`
for `++dist`'s PTRS branch, and on-ship-computed regression values plus
moment tests elsewhere, all on-ship (`-test %/tests/lib/rand ~`, all 91
green).

**Spec correction found during milestone 4**: `rand-spec.md` section 3.3
originally said PCG64 outputs the current state's permutation and THEN
advances. That's backwards -- pcg-c's actual reference
(`pcg_setseq_128_xsl_rr_64_random_r`, cross-checked against NumPy's vendored
copy) advances the state FIRST and outputs from the new state. Both the
spec and `++pcg` now follow the verified reference order; see the struck-
through note left in the spec for the record.

**Two real algorithm bugs caught during milestone 6's on-ship testing**
(both are exactly why this library tests against ship-computed values,
not hand-derivation): `++sample`'s alias-table `+build-loop` had the
small/large worklist assignment reversed (an index whose weight said it
belonged on one worklist was pushed onto the other); `++dist`'s
`+binomial` compared the uniform draw against the individual pmf term
instead of the accumulated CDF, which either returned 0 too often or
(for larger n) ran the recurrence past n and crashed with
subtract-underflow. Both are now documented in the arms' own doccords so
a future reader sees why the code looks the way it does.

Note on `++dist`: `+chi2`/`+student-t` aren't in milestone 5's literal arm
list but are one-line compositions of `+gamma`/`+normal` with no new
machinery, so they were added alongside rather than left in limbo until an
unscheduled milestone. `+categorical` has no `++rs` mirror: alias-table's
`.prob` is fixed at `@rd`, so it's precision-invariant.

Remaining milestones, in dependency order (see `rand-spec.md` section 13 for
full detail):

7. Non-float adapters (twoc, fixed, complex, posit lattice/unit).
8. Saloon `+rand-ray` (counter-window design) + replay-equality tests.
9. Full README rewrite + this file's final pass (ziggurat, BTPE, buffered
   Philox, posit rays, quire Monte Carlo note, `@rh`/`@rq` distributions —
   all deliberately deferred out of v1, per `rand-spec.md`).
