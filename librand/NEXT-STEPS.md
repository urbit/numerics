# `/lib/rand` — next steps

Status as of 2026-07-12. Milestones 1-4 (`rand-spec.md` section 13) are done:
`++split-mix`, `++philox`, `++seed`, `+step`, `+fork`, `++gen`, `++uni`,
`++pcg` -- KAT-verified against Vigna's reference SplitMix64, the Random123
`kat_vectors` file, an independent Python IEEE-754 encoder for `++uni`'s
float arms, and pcg-c's own seed=42/seq=54 demo convention for `++pcg`, all
on-ship (`-test %/tests/lib/rand ~`, all 40 green). `+step`'s `%pcg` branch
and `++uni`'s float/bits/below arms are now fully exercised on all three
engines (the earlier crash-stub for `%pcg` in `+step` is gone).

**Spec correction found during this milestone**: `rand-spec.md` section 3.3
originally said PCG64 outputs the current state's permutation and THEN
advances. That's backwards -- pcg-c's actual reference
(`pcg_setseq_128_xsl_rr_64_random_r`, cross-checked against NumPy's vendored
copy) advances the state FIRST and outputs from the new state. Both the
spec and `++pcg` now follow the verified reference order; see the struck-
through note left in the spec for the record.

Remaining milestones, in dependency order (see `rand-spec.md` section 13 for
full detail):

5. `++dist` at `@rd`: normal, expon, gamma, beta, bernoulli, geometric.
6. `++sample`: shuffle/permutation/choice/alias. `++dist` categorical/
   poisson/binomial/dirichlet, `@rs` variants.
7. Non-float adapters (twoc, fixed, complex, posit lattice/unit).
8. Saloon `+rand-ray` (counter-window design) + replay-equality tests.
9. Full README rewrite + this file's final pass (ziggurat, BTPE, buffered
   Philox, posit rays, quire Monte Carlo note, `@rh`/`@rq` distributions —
   all deliberately deferred out of v1, per `rand-spec.md`).
