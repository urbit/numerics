# `/lib/rand` and `/lib/i754rand` for Urbit

Deterministic pseudo-random number generation for numerical work (Monte
Carlo, ML init, shuffles, simulation). **Not cryptographic** — entropy
acquisition is Arvo's job (`eny`), cryptographic randomness is Zuse's job.

Full design is in `rand-spec.md` (repo root). Hoon reference implementation
first, jets follow (see `rand-spec.md` section 11).

## Two libraries, not one

`/lib/rand` is the **plumbing**: the three engines (`++philox`,
`++split-mix`, `++pcg`), seeding (`++seed`), the generic engine-dispatched
draw/fork (`+step`, `+fork`), the door facade (`++gen`), integer-only
uniform deviates (`++uni`: `+bits`/`+below`/`+between` — no floats), and
generic sampling (`++sample`: shuffle/permutation/choice/reservoir). The
`+$rng` type and its component engine-state types (`+$phil`/`+$sm64`/
`+$pcg64`) live in `/sur/rand`, not in the library itself, so any file
that only needs the type can get it via a lightweight `/-` import.

`/lib/i754rand` is the **porcelain**: IEEE-754 float generation
(`++uni`: `+rs`/`+rd`/`+rh`/`+rq`/`+rs-oo`/`+rd-oo`), Vose's alias method
(`++alias` — moved here from `/lib/rand`'s `++sample` because `+draw`
needs an actual uniform `@rd` float draw, so it isn't float-free the way
the rest of `++sample` is), and the distributions built on both
(`++dist`). This mirrors how `/lib/twoc`/`fixed`/`complex`/`unum` are
already kept separate from `/lib/math` in this codebase — `i754rand` is
`/lib/rand`'s own "math.hoon". Non-float output adapters (`twocrand`,
`fixedrand`, `complexrand`, `unumrand` — see `NEXT-STEPS.md`) are
`i754rand`'s siblings, not its dependents: none of them need floats for
their own core uniform-generation arms, only optionally for a "sample at
`@rd`, then quantize/convert" pattern (documented, not shipped as
dedicated wrapper arms — see `NEXT-STEPS.md`).

## Status (milestones 1-7 pass 1 of 9, `rand-spec.md` section 13)

Done, in `/lib/rand`:

- `++split-mix` — SplitMix64 (`+next`, `+split`, `+finalize`), KAT-checked
  against Vigna's reference C.
- `++philox` — Philox4x32-10 (`+block`, `+next`), the primary counter-based
  engine, KAT-checked against the three Random123 `kat_vectors` entries.
- `++seed` — `+from-atom`, `+from-eny`, `+fold-wide`, `+mix` (the two-word
  compression primitive `+fork` uses).
- `+step` — generic engine-dispatched draw across all three engines.
- `+fork` — path-sensitive key derivation across all three engine shapes.
- `++gen` — thin door facade wrapping `+step`/`+fork`.
- `++uni` — `+bits`, `+below` (Lemire, unbiased), `+between`.
- `++pcg` — PCG64 XSL-RR (`+next`, `+advance`, `+jump`), KAT-checked
  against pcg-c's own seed=42/seq=54 demo convention. **Corrects a spec
  bug found during implementation**: `rand-spec.md` originally said PCG
  outputs from the *current* state and advances after; the real reference
  (`pcg-c`, cross-checked against NumPy's vendored copy) advances the
  state FIRST and outputs from the new state. The spec and this
  implementation both now follow the verified reference order.
- `++sample` — `+shuffle`/`+permutation`/`+choice`/`+choices`/
  `+sample-n`/`+reservoir` (Algorithm R).

Done, in `/lib/i754rand`:

- `++uni` — the four float auras `+rs`/`+rd`/`+rh`/`+rq` plus open-open
  `+rs-oo`/`+rd-oo`, all exact bit constructions checked against an
  independent Python IEEE-754 encoder.
- `++alias` — Vose's alias method (`+build` + `+draw`) per rand-spec.md
  section 6.1, moved here from `++sample`.
- `++dist`, nested `++rd` (reference) / `++rs` (single precision, a
  mechanical re-instantiation of the same algorithms, per rand-spec.md's
  "each arm exists at @rd (reference) and @rs") — `+normal` (Marsaglia
  polar method), `+normal-mv`, `+expon` (inversion), `+gamma`
  (Marsaglia-Tsang, both alpha>=1 and the alpha<1 boost path), `+beta`,
  `+chi2`, `+student-t`, `+bernoulli`, `+geometric`, `+categorical` (rd
  only — thin wrapper over `++alias`'s table, which is fixed at @rd so
  has no meaningful rs variant), `+poisson` (Knuth for lambda<10,
  Hörmann's PTRS for lambda>=10 — verified against NumPy's
  `random_poisson_ptrs`), `+binomial` (inversion by CDF accumulation,
  crashes above `n*min(p,1-p) >= 30` where BTPE would be needed),
  `+dirichlet`. Moment tests (mean/variance regression at 50k draws,
  fixed seed) for normal/expon/gamma.

Done, non-float output adapters (rand-spec.md section 12):

- `/lib/twocrand` — `+twoc-full` (raw-bit passthrough), `+twoc-between`
  (Lemire-unbiased inclusive range in two's-complement order, via
  `/lib/twoc`'s width-keyed `+twid` door).
- `/lib/fixedrand` — `+fixed`, `+fixed-unit` (both raw-bit passthroughs —
  a fixed-point lattice is uniform by construction), `+fixed-between`
  (delegates to `twocrand`'s `+twoc-between`). The "sample at `@rd`,
  quantize" distribution pattern is documented, not shipped as dedicated
  wrapper arms; `tests/lib/fixedrand.hoon` proves the composition
  end to end via `/lib/fixed`'s `+from-rd` (added alongside this adapter).
- `/lib/complexrand` — `+cuniform`, `+normal-parts`, `+cnormal`,
  `+on-circle`, `+in-disk`, one arm set per component-width door
  (`+cd` reference precision, `+cs` mirror). Every arm draws floats
  directly, so this adapter depends on `/lib/i754rand` and `/lib/math`
  as well as `/lib/complex`. Building `+cnormal:cs` caught a real bug in
  `/lib/math`'s `@rs` `+invsqt2` constant (fixed, see `NEXT-STEPS.md`).

This completes milestone 7 pass 1. Not yet implemented: `/lib/unumrand`
(milestone 7 pass 2), the Saloon `+rand-ray` extension. See
`NEXT-STEPS.md` and `rand-spec.md` section 13 for the full milestone
order.

## Layout

```
librand/
  README.md
  NEXT-STEPS.md
  desk/sur/rand.hoon              :: +$rng, +$phil, +$sm64, +$pcg64
  desk/lib/rand.hoon              :: the plumbing
  desk/lib/i754rand.hoon          :: the IEEE-754 porcelain
  desk/lib/twocrand.hoon          :: two's-complement integer adapter
  desk/lib/fixedrand.hoon         :: fixed-point adapter
  desk/lib/complexrand.hoon       :: complex-number adapter
  desk/tests/lib/rand.hoon        :: -test %/tests/lib/rand ~
  desk/tests/lib/i754rand.hoon    :: -test %/tests/lib/i754rand ~
  desk/tests/lib/twocrand.hoon    :: -test %/tests/lib/twocrand ~
  desk/tests/lib/fixedrand.hoon   :: -test %/tests/lib/fixedrand ~
  desk/tests/lib/complexrand.hoon :: -test %/tests/lib/complexrand ~
```
