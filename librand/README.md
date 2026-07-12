# `/lib/rand` for Urbit

Deterministic pseudo-random number generation for numerical work (Monte
Carlo, ML init, shuffles, simulation). **Not cryptographic** — entropy
acquisition is Arvo's job (`eny`), cryptographic randomness is Zuse's job.

Full design is in `rand-spec.md` (repo root). Hoon reference implementation
first, jets follow (see `rand-spec.md` section 11).

## Status (milestones 1-6 of 9, `rand-spec.md` section 13)

Done:

- `++split-mix` — SplitMix64 (`+next`, `+split`, `+finalize`), KAT-checked
  against Vigna's reference C.
- `++philox` — Philox4x32-10 (`+block`, `+next`), the primary counter-based
  engine, KAT-checked against the three Random123 `kat_vectors` entries.
- `++seed` — `+from-atom`, `+from-eny`, `+fold-wide`, `+mix` (the two-word
  compression primitive `+fork` uses).
- `+step` — generic engine-dispatched draw across all three engines.
- `+fork` — path-sensitive key derivation across all three engine shapes.
- `++gen` — thin door facade wrapping `+step`/`+fork`.
- `++uni` — `+bits`, `+below` (Lemire, unbiased), `+between`, and the four
  float auras `+rs`/`+rd`/`+rh`/`+rq` plus open-open `+rs-oo`/`+rd-oo`, all
  exact bit constructions checked against an independent Python IEEE-754
  encoder.
- `++pcg` — PCG64 XSL-RR (`+next`, `+advance`, `+jump`), KAT-checked
  against pcg-c's own seed=42/seq=54 demo convention. **Corrects a spec
  bug found during implementation**: `rand-spec.md` originally said PCG
  outputs from the *current* state and advances after; the real reference
  (`pcg-c`, cross-checked against NumPy's vendored copy) advances the
  state FIRST and outputs from the new state. The spec and this
  implementation both now follow the verified reference order.

- `++dist`, now nested `++rd` (reference) / `++rs` (single precision, a
  mechanical re-instantiation of the same algorithms, per rand-spec.md's
  "each arm exists at @rd (reference) and @rs") — `+normal` (Marsaglia
  polar method), `+normal-mv`, `+expon` (inversion), `+gamma`
  (Marsaglia-Tsang, both alpha>=1 and the alpha<1 boost path), `+beta`,
  `+chi2`, `+student-t`, `+bernoulli`, `+geometric`, `+categorical` (rd
  only — thin wrapper over `++sample`'s alias table, which is fixed at
  @rd so has no meaningful rs variant), `+poisson` (Knuth for lambda<10,
  Hörmann's PTRS for lambda>=10 — verified against NumPy's
  `random_poisson_ptrs`), `+binomial` (inversion by CDF accumulation,
  crashes above `n*min(p,1-p) >= 30` where BTPE would be needed),
  `+dirichlet`. Moment tests (mean/variance regression at 50k draws,
  fixed seed) for normal/expon/gamma.
- `++sample` — `+shuffle`/`+permutation`/`+choice`/`+choices`/
  `+sample-n`/`+reservoir` (Algorithm R), and Vose's alias method
  (`++alias`: `+build` + `+draw`) per rand-spec.md section 6.1.

Not yet implemented: non-float adapters (twoc/fixed/complex/posit), the
Saloon `+rand-ray` extension. See `NEXT-STEPS.md` and `rand-spec.md`
section 13 for the full milestone order.

## Layout

```
librand/
  README.md
  NEXT-STEPS.md
  desk/lib/rand.hoon         :: the library
  desk/tests/lib/rand.hoon   :: -test %/tests/lib/rand ~
```
