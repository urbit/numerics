# `/lib/rand` — next steps

Status as of 2026-07-11. Milestone 1 (`rand-spec.md` section 13) is done:
`++split-mix` + `++seed`, KAT-verified against Vigna's reference SplitMix64
on-ship (`-test %/tests/lib/rand ~`, all green).

Remaining milestones, in dependency order (see `rand-spec.md` section 13 for
full detail):

2. `++philox` (Philox4x32-10, the primary counter-based engine) + Random123
   KAT vectors. `++fork` across all engines + path-sensitivity tests.
   `++gen` door facade.
3. `++uni` (bits/below/between/floats) + bias tests.
4. `++pcg` (PCG64 XSL-RR) + jump.
5. `++dist` at `@rd`: normal, expon, gamma, beta, bernoulli, geometric.
6. `++sample`: shuffle/permutation/choice/alias. `++dist` categorical/
   poisson/binomial/dirichlet, `@rs` variants.
7. Non-float adapters (twoc, fixed, complex, posit lattice/unit).
8. Saloon `+rand-ray` (counter-window design) + replay-equality tests.
9. Full README rewrite + this file's final pass (ziggurat, BTPE, buffered
   Philox, posit rays, quire Monte Carlo note, `@rh`/`@rq` distributions —
   all deliberately deferred out of v1, per `rand-spec.md`).
