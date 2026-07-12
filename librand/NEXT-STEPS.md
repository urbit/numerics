# `/lib/rand` — next steps

Status as of 2026-07-11. Milestones 1-2 (`rand-spec.md` section 13) are done:
`++split-mix`, `++philox`, `++seed`, `+step`, `+fork`, `++gen` -- KAT-verified
against Vigna's reference SplitMix64 and the Random123 `kat_vectors` file
on-ship (`-test %/tests/lib/rand ~`, all 23 green).

Note for whoever picks up milestone 4: `+step`'s `%pcg` branch currently
crashes (`~|  %rand-pcg-step-not-yet-implemented`) since a real draw needs
PCG's xsl-rr output permutation, which doesn't exist yet. `+fork`'s `%pcg`
branch is already fully implemented (forking only remixes state/inc via
`+mix`, no draw logic needed), so only `+step` needs a fix once `++pcg` lands.

Remaining milestones, in dependency order (see `rand-spec.md` section 13 for
full detail):

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
