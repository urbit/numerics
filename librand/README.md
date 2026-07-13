# `/lib/rand`: deterministic pseudo-random number generation for Urbit

Reproducible random number generation for numerical work — Monte Carlo,
ML initialization, shuffles, simulation. **Not cryptographic**: entropy
acquisition is Arvo's job (`eny`), cryptographic randomness is Zuse's job.
Given the same seed, every arm here produces the same output on any
ship, any time, forever — that's the entire point.

Full design is in `rand-spec.md` (repo root). All nine of its milestones
are done: three engines, integer and float uniform deviates, a full
distribution suite, sampling utilities, four non-float output adapters,
and a Lagoon array-filling layer in Saloon. Hoon reference implementation
first; jets are a follow-up (`rand-spec.md` section 11).

## Layout: seven libraries, not one

The whole surface is split by concern, mirroring how `/lib/twoc`/`fixed`/
`complex`/`unum` are already kept separate from `/lib/math` in this
codebase:

- **`/sur/rand`** — just the types (`+$rng`, `+$phil`, `+$sm64`,
  `+$pcg64`), so anything that only needs the type gets it via a
  lightweight `/-` import instead of pulling in the whole library.
- **`/lib/rand`** — the **plumbing**: the three engines, seeding,
  engine-dispatched draw/fork, integer-only uniform deviates, generic
  sampling. Every other library here imports this one.
- **`/lib/i754rand`** — the **porcelain**: IEEE-754 float generation,
  Vose's alias method, and the full distribution suite. This is
  `/lib/rand`'s own "math.hoon" — kept separate so the plumbing and the
  non-float adapters below never need to pull in `/lib/math`.
- **`/lib/twocrand`**, **`/lib/fixedrand`**, **`/lib/complexrand`**,
  **`/lib/unumrand`** — non-float output adapters (two's-complement
  integers, fixed-point, complex numbers, posits). Siblings of
  `i754rand`, not dependents of it: only `complexrand` needs floats for
  its own uniform-generation arms (every arm draws floats directly); the
  others only optionally reach for `i754rand` to compose a "sample at
  `@rd`, then quantize/convert" distribution (documented as a pattern,
  not shipped as dedicated wrapper arms — see `NEXT-STEPS.md`).
- **`/lib/saloon`** (a different desk) — the array-filling layer,
  `+rand-ray`, is not here; it lives directly in Saloon's `+sa` core (see
  below).

## What's here

### `/lib/rand` (plumbing)

- **Engines**: `++philox` (Philox4x32-10, counter-based, the primary
  engine — its whole design point is that element `i` of an array only
  needs counter `ctr0+i`, so a jet can fill array elements in parallel
  and land bit-identical results regardless of thread scheduling),
  `++split-mix` (SplitMix64, a small sequential generator and the
  seed-mixing primitive every engine's seeding/`+fork` path depends on),
  `++pcg` (PCG64 XSL-RR, with `O(log n)` `+jump`). All three KAT-checked
  against independent reference implementations (Random123's
  `kat_vectors`, Vigna's public-domain SplitMix64 C, pcg-c's own
  seed=42/seq=54 demo convention).
- **`++seed`** — `+from-atom`, `+from-eny`, `+fold-wide`, `+mix`.
- **`+step`/`+fork`** — the generic engine-dispatched draw and the
  path-sensitive key derivation (JAX-style: fork by path instead of
  threading one sequential state, so a draw inserted in one branch of a
  computation never perturbs a sibling branch) across all three engine
  shapes. `++gen` is a thin door facade over both.
- **`++uni`** (integer only) — `+bits`, `+below` (Lemire's unbiased
  method — the primitive every bounded-integer draw in this codebase
  uses; modulo bias is never an option here), `+between`.
- **`++sample`** — `+shuffle`/`+permutation`/`+choice`/`+choices`/
  `+sample-n`/`+reservoir` (Algorithm R).

### `/lib/i754rand` (porcelain)

- **`++uni`** — the four float auras (`+rs`/`+rd`/`+rh`/`+rq`) plus
  open-open variants (`+rs-oo`/`+rd-oo`), exact bit constructions checked
  against an independent Python IEEE-754 encoder.
- **`++alias`** — Vose's alias method (`+build`/`+draw`), moved here
  from `++sample` because `+draw` needs an actual float draw to decide
  accept-vs-redirect.
- **`++dist`**, nested `++rd` (reference, double precision) / `++rs`
  (single precision, a mechanical re-instantiation of the same
  algorithms) — `+normal` (Marsaglia polar), `+normal-mv`, `+expon`
  (inversion), `+gamma` (Marsaglia-Tsang), `+beta`, `+chi2`,
  `+student-t`, `+bernoulli`, `+geometric`, `+categorical` (`@rd`
  only — the alias table's probabilities are fixed at that precision),
  `+poisson` (Knuth for `lambda<10`, Hörmann's PTRS above — verified
  against NumPy's `random_poisson_ptrs`), `+binomial` (inversion by CDF
  accumulation; crashes above `n*min(p,1-p) >= 30`, where BTPE would be
  needed — see `NEXT-STEPS.md`), `+dirichlet`. Moment-tested (mean/
  variance regression at 50k draws, fixed seed) for normal/expon/gamma.
  `@rh`/`@rq` distributions are out of scope for v1 (see `NEXT-STEPS.md`).

### The four non-float adapters (`rand-spec.md` section 12)

- **`/lib/twocrand`** — `+twoc-full` (raw-bit passthrough) and
  `+twoc-between` (Lemire-unbiased inclusive range in two's-complement
  order, via `/lib/twoc`'s width-keyed `+twid` door). Depends only on
  `/lib/rand` + `/lib/twoc`.
- **`/lib/fixedrand`** — `+fixed`, `+fixed-unit` (both raw-bit
  passthroughs — a fixed-point lattice is uniform by construction),
  `+fixed-between` (delegates to `twocrand`). Needed a new `+from-rd`
  added to `/lib/fixed`, mirroring its existing `+from-rs`.
- **`/lib/complexrand`** — `+cuniform`, `+normal-parts`, `+cnormal`
  (deliberately distinct from `+normal-parts`: "complex Gaussian" means
  different things to signal-processing and statistics users),
  `+on-circle`, `+in-disk`. One arm set per component-width door (`+cd`
  reference precision, `+cs` mirror), matching `/lib/complex`'s own ship
  order. Every arm draws floats directly, so this is the one adapter
  that depends on `/lib/i754rand`, `/lib/math`, and `/lib/complex`.
- **`/lib/unumrand`** — two *inequivalent* uniform semantics over posits,
  named so they can't be confused (posits are tapered: consecutive bit
  patterns are not evenly spaced in value). `+posit-lattice` is uniform
  over raw bit patterns excluding NaR (a fuzzing/property-testing
  primitive — the induced value distribution is only approximately
  log-uniform); it ships at all five width doors (`rpb`/`rph`/`rps`/
  `rpd`/`rpq`), with no bit-count subtlety at any width. `+posit-unit` is
  uniform over *values* on `[0,1)`, exact; it's scoped to posit8/16/32
  only (see `NEXT-STEPS.md` for the `k=4n` bit-count derivation behind
  its exactness claim, and why posit64/128 aren't included). Needs no new
  `/lib/unum` plumbing — `+from-rh/rs/rd/rq` already exist at every width
  door. Verified against a new exact-rational oracle
  (`tools/posit_unit_check.py`) via chi-square at posit8 (100,000
  ship-drawn draws, p=0.108) plus bit-exact cross-checks at every
  in-scope width.

None of the four ship dedicated per-distribution wrapper arms for the
"sample at `@rd`, quantize/convert" pattern (fixed-point and posit
distributions) — that would be a dozen-plus nearly-identical one-liners
composing two already-existing arms. Instead each adapter documents the
one-line composition and has one test proving it works end to end.

### `+rand-ray` (in `/lib/saloon`, `rand-spec.md` section 8)

Not part of this library's own file tree — it's a new `+|  %rand` section
in Saloon's `+sa` core, filling a Lagoon `$ray` element-by-element in
row-major (C) order: `+fill-uniform`, `+fill-normal`, `+fill-expon`
(`%i754` rays, bloq 5/6 — `@rs`/`@rd` only; posit rays are deferred),
`+fill-below` (`%uint` rays via Lemire). Philox elements get per-element
counter treatment — the entire reason Philox is the primary engine:
`+fill-uniform`'s single non-rejecting draw per element assigns element
`i` counter `ctr0+i` directly, so the whole fill decomposes into `n`
independent draws a future jet can parallelize. The rejection-based fills
(`+fill-normal`/`+fill-expon`/`+fill-below`) get a wider per-element
counter *window* instead (`ctr0 + i*2^32`): each element's rejection loop
walks freely within its own window, and the returned rng's counter is
forced to `ctr0 + n*2^32` regardless of how many sub-draws any element
actually used, so the post-state is a pure function of `n` — with an
explicit crash if a single element's rejection loop ever walks past its
own window (astronomically improbable, but checked rather than silently
overflowing into the next element's window). Non-Philox engines have no
comparable jump primitive and aren't the ones a jet would parallelize
anyway, so they just thread sequentially.

Every arm is bit-exact cross-checked against a direct call to the
underlying `/lib/rand`/`/lib/i754rand` primitive at the expected counter
— see `saloon/desk/tests/lib/saloon-rand-ray.hoon` (in the `saloon` desk,
not here).

## Testing philosophy

Every arm is verified against ship-computed output before being trusted —
this codebase does not hand-derive expected values for anything involving
floating-point transcendentals, bit-level float/posit construction, or
RNG state transitions, because a hand or Python derivation isn't
guaranteed to match this codebase's own kernels bit-for-bit. Known-answer
tests come from independent reference implementations wherever one
exists (Random123, Vigna's SplitMix64, pcg-c, NumPy's PTRS, an exact-
rational posit oracle); everything else is a ship-computed regression
value pinned into the test file. Two real algorithm bugs (a reversed
worklist in the alias method, a wrong-accumulator comparison in
`+binomial`) and two real pre-existing bugs in dependencies (a malformed
`@rs` float literal in `/lib/math`'s `+invsqt2`; a sample-shape mismatch
in Saloon's scalar transcendental dispatch) were caught exactly this way.

```
-test %/tests/lib/rand ~
-test %/tests/lib/i754rand ~
-test %/tests/lib/twocrand ~
-test %/tests/lib/fixedrand ~
-test %/tests/lib/complexrand ~
-test %/tests/lib/unumrand ~
-test %/tests/lib/saloon-rand-ray ~     :: in the saloon desk
```

## Layout

```
librand/
  README.md
  NEXT-STEPS.md
  rand-spec.md                    :: (repo root) the governing design spec
  desk/sur/rand.hoon               :: +$rng, +$phil, +$sm64, +$pcg64
  desk/lib/rand.hoon                :: the plumbing
  desk/lib/i754rand.hoon            :: the IEEE-754 porcelain
  desk/lib/twocrand.hoon            :: two's-complement integer adapter
  desk/lib/fixedrand.hoon           :: fixed-point adapter
  desk/lib/complexrand.hoon         :: complex-number adapter
  desk/lib/unumrand.hoon            :: posit adapter
  desk/tests/lib/rand.hoon          :: -test %/tests/lib/rand ~
  desk/tests/lib/i754rand.hoon      :: -test %/tests/lib/i754rand ~
  desk/tests/lib/twocrand.hoon      :: -test %/tests/lib/twocrand ~
  desk/tests/lib/fixedrand.hoon     :: -test %/tests/lib/fixedrand ~
  desk/tests/lib/complexrand.hoon   :: -test %/tests/lib/complexrand ~
  desk/tests/lib/unumrand.hoon      :: -test %/tests/lib/unumrand ~
  tools/posit_unit_check.py         :: +posit-unit's exact-rational oracle

saloon/desk/lib/saloon.hoon         :: +rand-ray lives in +sa's +| %rand section
saloon/desk/tests/lib/saloon-rand-ray.hoon
```

See `NEXT-STEPS.md` for what's deliberately deferred past v1 (ziggurat,
BTPE, buffered Philox, posit rays, quire Monte Carlo, `@rh`/`@rq`
distributions) and the full development history/decision log.
