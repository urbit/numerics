# SPEC: `/lib/rand` — Deterministic Random Number Generation for Urbit Numerics

**Target repo:** `urbit/numerics`, new directory `librand/` (structured like `libmath/`), with a Saloon extension for ray-filling.
**Status:** Draft for implementation. Hoon reference implementation first; jets follow.
**Dependencies:** `/lib/math` (float transcendentals), `/lib/twoc` (width-keyed modular integers), optionally `/lib/unum` and `/lib/fixed` for output adapters. Saloon layer depends on `/lib/lagoon`.

---

## 0. Design principles (read before writing code)

1. **The Hoon is the spec.** Every generator and every distribution transform must be defined such that a C jet can reproduce it **bit-for-bit**. This is the same discipline as `/lib/math`'s SoftFloat parity. Consequences:
   - All float arithmetic in distribution transforms goes through the stdlib `@r*` doors or `/lib/math` kernels (which have jet-parity guarantees). No algorithm may depend on extended-precision intermediates, FMA availability, or evaluation order beyond what the Hoon states.
   - All integer arithmetic is exact `@` arithmetic masked to width via `/lib/twoc` (`+twid` at width 32/64/128) or explicit `(mod . (bex w))`. Never rely on implicit truncation.

2. **Counter-based first.** The primary generator is **Philox4x32-10** (Salmon et al., "Parallel Random Numbers: As Easy as 1, 2, 3", SC'11; the Random123 library is the reference). Rationale:
   - Output is a pure function `(key, counter) -> 4 x u32`. No sequential state dependency, so a Lagoon jet can fill a ray in parallel and still match the sequential Hoon loop exactly.
   - Reproducibility across ships, runtimes, and replays (Nockchain verification contexts) falls out for free.
   - Stream splitting is trivial and collision-free (distinct keys = distinct streams).

3. **Sequential generators are provided but secondary.** `+split-mix` (SplitMix64) for seed expansion and cheap draws; `+pcg` (PCG64 XSL-RR 128/64) for users who want a conventional stateful generator. **xoshiro is explicitly excluded** — nothing it offers matters here and its low-bit weaknesses are a support liability.

4. **Not cryptographic. Say so loudly.** This library produces *deterministic pseudo-randomness for numerical work* (Monte Carlo, ML init, shuffles, simulation). It must NOT be used for key material, nonces, or anything adversarial. Entropy acquisition is Arvo's job (`eny`); cryptographic randomness is Zuse's job. The library docstring, the README, and each seeding arm that accepts `eny` must carry this warning.

5. **State threading is explicit and total.** Every drawing arm returns `[value new-state]` over plain-data state (section 2.1 defines the full discipline: functional core, door facade, key derivation). No arm hides state mutation. Rejection loops are fine (they thread state a variable number of times) but must be documented as variable-consumption arms. Relationship to stdlib `+og`: we adopt its door *ergonomics* as a facade, not its engine (iterated SHA-256, slow, unjetted for this purpose) and not its historical misuse (persisting the core).

6. **Match house style.** Doors and arm docs follow the existing repo conventions (see `/lib/math`): `::  +arm:  type -> type` header comment, `Examples` block with dojo output, `Source` marker, `~/` jet hints with `~%  %non  ..part  ~` registration at the core.

---

## 1. Module layout

```
librand/
  README.md
  NEXT-STEPS.md
  desk/lib/rand.hoon        :: everything below except the Saloon layer
saloon/desk/lib/saloon.hoon :: gains a +rand-ray core (section 7)
```

Single library file, one top-level core, sub-cores per concern:

```
/lib/rand
|%
++  philox    :: counter-based engine (primary)
++  split-mix :: SplitMix64: seed expansion + cheap sequential generator
++  pcg       :: PCG64 (XSL-RR 128/64) sequential generator
++  seed      :: seeding utilities (from @, from eny, hash-mixing)
++  uni       :: uniform deviates: bits, bounded ints, floats per aura
++  dist      :: nonuniform distributions (normal, exp, gamma, ...)
++  sample    :: shuffles, permutations, choice, reservoir, alias tables
--
```

---

## 2. Core types

```hoon
::  Philox key/counter state.  ctr is the 128-bit counter as a single @,
::  key is the 64-bit key as a single @.  Both are plain atoms; width
::  discipline is enforced by masking, never by aura tricks.
+$  phil  [key=@ ctr=@]

::  Sequential generator state.
+$  sm64  @        :: SplitMix64: 64-bit state
+$  pcg64  [state=@ inc=@]   :: 128-bit state, 128-bit odd increment

::  A generic stream: tagged union so distribution code is engine-agnostic.
+$  rng
  $%  [%phil p=phil]
      [%sm64 s=sm64]
      [%pcg p=pcg64]
  ==
```

Distribution and sampling arms take and return `rng`, calling a single internal `+step` gate, which lives at the top level of `/lib/rand` (not nested under any one engine core, since it dispatches across the `+$rng` tagged union):

```hoon
::  +step: rng -> [u64=@ rng]   (one 64-bit draw, engine-dispatched)
```

This keeps `+dist` and `+sample` written once, not per-engine. Resolved design (no buffering, no caching — simplicity and an unambiguous spec beat throughput at v1): for `%phil`, `+step` calls `next:philox`, which returns the full 128-bit block `out = (block key ctr)`; `+step` takes **bits `[0,64)` of that atom**, i.e. lanes `c0'` and `c1'` (the low two lanes of the little-endian repack — `out = c3'*2^96 + c2'*2^64 + c1'*2^32 + c0'`, so bits `[0,64)` = `c1'*2^32 + c0'`), and increments the counter by 1. The top 64 bits (`c2'`, `c3'`) are discarded; wasting them is acceptable at v1, and a buffered variant that reuses both halves is a documented follow-up in NEXT-STEPS. For `%sm64`, `+step` calls `next:split-mix` directly (already 64 bits). For `%pcg`, `+step` calls `next:pcg`, whose `xsl-rr` output is already a 64-bit rotr64 result.

### 2.1 State discipline and API surface

Three sanctioned usage patterns; the API serves all three.

**(a) Functional core — canonical.** Every drawing arm is `[value new-rng]` over the plain-data `+$rng` noun. Call-site idiom is `=^`:

```hoon
=^  x    rng  (rd:uni:rand rng)
=^  y    rng  (normal:dist:rand rng)
```

Every arm's `Examples` block demonstrates the `=^` form. Agent state stores the `+$rng` **noun**.

**(b) Door facade — ergonomic sugar, never a storage format.** Provide `++gen`, an `+og`-style door whose sample is the `rng` noun and whose arms mirror the functional core, returning `[value _+>.$]`:

```hoon
=/  g  ~(. gen:rand (from-atom:seed:rand %phil 42))
=^  x  g  rd:g
```

**HARD RULE, stated in the door's docstring:** the door is a call-site convenience. Persisting it in agent state pins a battery across library upgrades — the classic `og` misuse. Persist the `+$rng` data noun; reconstruct the door locally. The facade adds no capability; it is a strict wrapper over (a), so jets register on the functional arms only.

**(c) Key derivation — JAX-style, first-class, the payoff of counter-based-first.**

```hoon
::  +fork:  [rng @] -> rng   (deterministic, independent child stream)
++  fork  |=([r=rng salt=@] ^-(rng))
```

Per engine: `%phil` — child key = `(mix key salt)`, the two-word compression pinned down in section 4, counter reset to 0. `%pcg` — child increment derived the same way (`(mix inc salt)`), forced odd; state re-mixed. `%sm64` — the SplitMix split construction. Distinct salts give computationally independent streams.

Instead of threading one sequential state through a computation tree, callers fork by path: per-layer keys for Maroon weight initialization, per-event keys in agents (`(fork base eny-of-event-id)`), per-cell keys in simulations. Properties sequential threading cannot provide: draws are independent of sibling evaluation order, and inserting a draw in one branch does not perturb any other branch. Document this as the *recommended* pattern for tree-structured and ML workloads; sequential threading (a) remains correct and simpler for linear code.

Nesting rule to spec: `(fork (fork r a) b)` must differ from `(fork (fork r b) a)` and from `(fork r (cat 6 a b))` — i.e., the mix is genuinely path-sensitive. Test this.


---

## 3. Engines

### 3.1 Philox4x32-10 (`++philox`)

Reference: Random123. Constants (u32):

```
M0 = 0xD2511F53   M1 = 0xCD9E8D57
W0 = 0x9E3779B9   W1 = 0xBB67AE85   (bumpkeys per round)
Rounds = 10
```

Arms:

- `++block  |=([key=@ ctr=@] ^-(@))` — the pure Philox function. Split `ctr` into `c0 c1 c2 c3` (u32 little-endian lanes), `key` into `k0 k1`. Ten rounds of the Philox S-box:
  ```
  hi0.lo0 = mulhilo32(M0, c0);  hi1.lo1 = mulhilo32(M1, c2)
  c0' = hi1 ^ c1 ^ k0 ;  c1' = lo1
  c2' = hi0 ^ c3 ^ k1 ;  c3' = lo0
  k0 += W0 ; k1 += W1   (mod 2^32)
  ```
  Return the four output lanes repacked as one 128-bit atom, same lane order.
  `mulhilo32` is exact: `p = (mul a b)`, `hi = (rsh [0 32] p)`, `lo = (end [0 32] p)`.
- `++next  |=(p=phil ^-([out=@ p=phil]))` — `out = (block key ctr)`, new state `[key (mod +(ctr) (bex 128))]`.
- **Test vectors:** the spec MUST embed the Random123 known-answer test: `philox4x32-10(ctr = {0,0,0,0}, key = {0,0})`, `ctr = key = all 0xffffffff...`, and the pi-digits vector from the Random123 `kat_vectors` file. These go in `/tests/lib/rand.hoon` and are non-negotiable — an engine that fails KAT is wrong, full stop.

### 3.2 SplitMix64 (`++split-mix`)

Reference: Steele, Lea, Flood (OOPSLA'14); Vigna's public-domain C.

```
next(s):  s += 0x9E3779B97F4A7C15                (mod 2^64)
          z = s
          z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9   (mod 2^64)
          z = (z ^ (z >> 27)) * 0x94D049BB133111EB   (mod 2^64)
          return z ^ (z >> 31)
```

Arms: `++next |=(s=sm64 [out=@ s=sm64])`, plus `++split` returning two decorrelated child states (gamma-stepping per the paper is overkill; two `+next` draws as child seeds is adequate and simpler — document the choice).

KAT: first 5 outputs from seed `0` and seed `0xDEADBEEF` against the reference C.

### 3.3 PCG64 XSL-RR 128/64 (`++pcg`)

Reference: O'Neill 2014, pcg-random.org. LCG multiplier `0x2360ED051FC65DA44385DF649FCCF645` (128-bit), user-chosen odd increment (stream id). Output: `xsl-rr`: `rot = state >> 122`; `xored = (state >> 64) ^ (state & mask64)`; output = rotr64(xored, rot).

**Correction (implementation phase, verified against `pcg-c`'s `include/pcg_variants.h` — `pcg_setseq_128_xsl_rr_64_random_r` — and independently against NumPy's vendored `pcg64.orig.h`, both consistent): an earlier draft of this spec had the operation order backwards.** The real reference does `pcg_setseq_128_step_r(rng); return pcg_output_xsl_rr_128_64(rng->state);` — i.e. **advance the state FIRST, then compute the output permutation from the NEW (already-advanced) state.** (The draft text below, retained struck through for the record, claimed the opposite: "State update AFTER output extraction (i.e., output the *current* state's permutation, then advance)" — that is wrong; do not implement it that way.) Per this spec's own closing rule (section 13): the reference implementation wins over the spec text. This is exactly the classic off-by-one-draw divergence the original text worried about, just resolved in the other direction.

Include `++jump` (advance by `2^64` steps via LCG skip-ahead: standard `O(log n)` modular matrix trick on `(a, c)`) for stream partitioning.

KAT: `pcg64_random_r` demo output from the reference distribution, seed `42, 54`.

---

## 4. Seeding (`++seed`)

- `++from-atom  |=(a=@ ^-(rng))` — hash `a` through SplitMix64 to fill whichever engine (parameterize: `|=([eng=?(%phil %sm64 %pcg) a=@] ...)`). Philox: key = first SplitMix output, ctr = 0. PCG: state/inc from two outputs, inc forced odd via `(con inc 1)`.
- `++from-eny  |=(eny=@uvJ ^-(rng))` — same pipeline, over Arvo entropy. `eny` is 512 bits; fold it to 64 bits first by chunking into eight 64-bit words (low to high) and folding sequentially through `+mix` (`acc = 0; acc = (mix acc chunk)` for each chunk in order) before handing the result to `+from-atom`'s pipeline. Carries the not-for-crypto warning.
- `++mix  |=([a=@ b=@] @)` — combine two 64-bit seed words (e.g., ship + per-agent salt) into one. **Fully specified as a two-word compression, not a single finalizer call on the raw 128-bit concatenation** (the finalizer below is a function of one 64-bit word; a 128-bit input must be folded, not fed in directly):
  ```
  w0 = a mod 2^64  ;  w1 = b mod 2^64
  z0 = finalize(w0)              :: the z-mixing tail of +next (3.2): z=z0; z=(z^(z>>30))*C1; z=(z^(z>>27))*C2; return z^(z>>31) — no `s += golden` step, that increment is +next's state-advance, not part of the finalizer
  return finalize(z0 XOR w1)
  ```
  Order-sensitive by construction: `(mix a b)` differs from `(mix b a)` in general, since `a` is absorbed through one extra `finalize` application relative to `b`. This is what makes `(fork (fork r a) b)` differ from `(fork (fork r b) a)` (section 2.1c) — each nested `fork` call adds another `finalize` layer at a different depth. Deterministic, no rng threading (pure `@ -> @`). For inputs wider than 64 bits (e.g. `eny`), chunk into 64-bit words and fold sequentially as shown above for `+from-eny` — do not truncate.
  **Naming note (implementation footgun):** this arm's name collides with the Hoon standard library's `++mix` (bitwise XOR). Inside `++seed`, an unqualified `(mix a b)` call resolves to this shadowing local arm; any internal use of stdlib XOR inside `++seed` must reference it explicitly (e.g. `^mix`). Callers outside `++seed` are unaffected. Note this in the arm's doc comment so nobody is bitten by it during implementation.

No implicit global seeding. There is no "default rng"; callers own their state.

---

## 5. Uniform deviates (`++uni`)

All arms `|=([r=rng ...] [out new-rng])` unless noted.

### 5.1 Integers

- `++bits  |=([r=rng n=@] [@ rng])` — `n` uniform bits, `n <= 64` per draw, looping for larger `n` (assemble little-endian, low draw first — specify this).
- `++below  |=([r=rng n=@] [@ rng])` — uniform in `[0, n)`, **unbiased**, via Lemire's multiply-shift with rejection (Lemire 2019, "Fast Random Integer Generation in an Interval"):
  ```
  draw x: u64;  m = x * n (exact @ arithmetic, 128-bit product)
  l = m mod 2^64
  if l < n:  t = (2^64 - n) mod n;  while l < t: redraw x, recompute
  return m >> 64
  ```
  Crash (`?>`) on `n = 0`. This is the primitive everything else (shuffle, choice) uses. Do NOT ship modulo-bias `(mod x n)` anywhere, including tests.
- `++between  |=([r=rng a=@s b=@s] [@s rng])` — inclusive signed range via `+below` on the width, offset. Uses `si` arithmetic; crash on `a > b`.

### 5.2 Floats

One arm per aura, matching the repo's four-precision pattern:

- `++rs  |=(r=rng [@rs rng])` — uniform in `[0,1)`: draw 24 bits, `out = (mul:rs (sun:rs bits) 2^-24)` — implemented as the exact bit construction `(bits << 0) * 0x1p-24`, i.e. build the float by `~(sun rs %n)` then multiply by the constant `0x3380.0000` (`2^-24`). Every representable output is a multiple of `2^-24`; distribution is exactly uniform over that lattice. Document that `0` is possible and `1` is not.
- `++rd` — 53 bits, times `2^-53` (`0x3CA0.0000.0000.0000`).
- `++rh` — 11 bits, times `2^-11`.
- `++rq` — 113 bits, times `2^-113`.
- `++rs-oo` / `++rd-oo` (open-open `(0,1)`): same but `(bits + 0.5) * 2^-w` via drawing `w` bits (the SAME bit count as the closed-open variant, not `w-1` — drawing `w-1` bits confines the output to `(0, 0.5)`, not `(0,1)`; this was an error in an earlier draft) and setting the low half-ulp: implement as `((2*bits + 1) * 2^-(w+1))` with exact integer construction. With `bits` ranging over `[0, 2^w)`, `2*bits+1` ranges over the odd integers in `[1, 2^(w+1))`, so the result spans `(0,1)` on a lattice of spacing `2^-w`, symmetric about 0.5, never touching either endpoint. Needed by Box-Muller/log-based transforms so `log(0)` never fires.
- Posit, fixed, twoc, and complex outputs: section 12. Posit spacing is non-uniform, so posit uniformity is a distinct design with two semantics — resolved there, not deferred.

---

## 6. Distributions (`++dist`)

All doors keyed the way `/lib/math` doors are where a rounding mode matters; internally force `%n` like the math kernels do and say so. Each arm exists at `@rd` (reference) and `@rs` (routed through the same algorithm at single precision). `@rh`/`@rq` deferred to NEXT-STEPS.

Algorithm selections — chosen for jet-parity determinism (comparison-and-arithmetic only, or transcendentals that route through `/lib/math`'s bit-specified kernels):

| Arm | Algorithm | Reference | Notes |
|---|---|---|---|
| `++normal` | **Polar Marsaglia** (Box-Muller polar variant) | Marsaglia & Bray 1964 | Rejection on `s >= 1` or `s == 0`; uses `+log`/`+sqt` from `/lib/math`. Returns ONE deviate, discards the pair-mate at v1 (caching the mate makes state opaque; document the waste). Ziggurat is a NEXT-STEPS optimization — it needs embedded tables that must be spec'd exactly, don't let Sonnet improvise them at v1. |
| `++normal-mv` | mean/sigma wrapper | — | `mu + sigma * z`. Crash on `sigma < 0`. |
| `++expon` | inversion: `-log(u)` over `(0,1)` draw | — | rate parameter `lambda`: divide. Crash `lambda <= 0`. |
| `++gamma` | **Marsaglia–Tsang** squeeze | Marsaglia & Tsang 2000 | `alpha >= 1` direct; `alpha < 1` via boost `gamma(alpha+1) * u^(1/alpha)` (uses `+pow` from math). Crash `alpha <= 0`. |
| `++beta` | two gammas: `x/(x+y)` | — | |
| `++chi2` | `gamma(k/2, 2)` | — | |
| `++student-t` | normal / sqrt(chi2/k) | — | |
| `++poisson` | Knuth product for `lambda < 10`; **PTRS** (Hörmann 1993) transformed rejection for `lambda >= 10` | Hörmann, "The transformed rejection method for generating Poisson random variables" | Knuth-only is O(lambda) and will bite someone; spec both branches and the exact switch point. |
| `++binomial` | inversion for `n*min(p,1-p) < 30`; **BTPE** deferred | Kachitvichyanukul & Schmeiser 1988 | v1 ships inversion + a documented crash (`?>`) above the threshold rather than a slow or wrong large-n path. BTPE in NEXT-STEPS. |
| `++geometric` | `ceil(log(u)/log(1-p))` | — | edge: `p = 1` returns 1; crash `p <= 0` or `p > 1`. |
| `++bernoulli` | `u < p` | — | returns `?`. |
| `++dirichlet` | k gammas, normalized | — | list in, list out. |
| `++categorical` | via `++alias` (section 6.1) | Vose 1991 | |

### 6.1 Alias method (`++alias` in `++sample`)

Vose's O(n) construction, O(1) draw. Table type:

```hoon
+$  alias-table  [n=@ prob=(list @rd) alias=(list @ud)]
```

Construction is pure (no rng); drawing takes `[t=alias-table r=rng]`. Weights in as `(list @rd)`, non-negative, at least one positive (crash otherwise). Normalization inside construction. Spec the small/large worklist algorithm precisely — it's the classic place for off-by-one and float-compare bugs; require the construction to use `@rd` `%n` arithmetic only.

---

## 7. Sampling utilities (`++sample`)

- `++shuffle  |=([r=rng l=(list)] [(list) rng])` — **Fisher–Yates**, iterating from the END down (`i` from `n-1` to `1`, swap with `(below r +(i))`). Implement over a `(map @ud *)` or flopped-list index structure to avoid O(n²) `snag`/`oust`; O(n log n) via map is fine, note the tradeoff. Wet gate (`|*`) so element type is preserved.
- `++permutation  |=([r=rng n=@] [(list @ud) rng])` — shuffle of `(gulf 0 (dec n))`.
- `++choice  |=([r=rng l=(list)] [* rng])` — one element, uniform. Crash on empty.
- `++choices` — n with replacement.
- `++sample-n` — n WITHOUT replacement: partial Fisher–Yates (first n of a permutation), not repeated rejection.
- `++reservoir  |=([r=rng n=@ l=(list)] ...)` — Algorithm R (Vitter). Included because agents streaming events want it.

---

## 8. Saloon layer (`+rand-ray`)

In `/lib/saloon`, a core taking a Lagoon `meta` and an `rng`, returning `[ray rng]`:

- `++fill-uniform  |=([=meta r=rng] [ray rng])` — uniform `[0,1)` at the meta's float aura/bloq. Iterates the counter per element in **row-major (C) order** — the order must be specified because the jet will parallelize and must land identical bits per element. For `%phil`, element `i` uses counter `ctr0 + i` and the post-state is `ctr0 + n`; this is the jet-parallelism payoff and the reason Philox is primary. For sequential engines, elements are drawn in row-major sequence.
- `++fill-normal`, `++fill-expon` — same pattern. **Constraint:** rejection-based transforms (polar normal) consume variable draws, which breaks per-element counter assignment. Resolve: per-element sub-counter space — element `i` owns counters `ctr0 + i*2^32 ..`, rejection walks within its window (window exhaustion is astronomically improbable; crash if it happens). Post-state after filling `n` elements is `ctr0 + n*2^32` — i.e. the base counter advances past the *entire* window block for every element, not just the sub-counters actually consumed within each element's rejection loop. This keeps the post-state a pure function of `n`, independent of how many draws each element's rejection happened to need, so replay and post-state composition (e.g. chaining another `+fill-*` call) stay simple. Spec this exactly; it is the one genuinely novel design element in the library.
- Integer rays: `++fill-below` for `%uint` rays via Lemire.

Aura/type support at v1: `%r32`/`%r64` rays. Posit rays deferred.

---

## 9. Errors and edge cases (uniform policy)

- Domain violations crash via `?>` with `~|` tags (`%rand-empty-list`, `%rand-bad-prob`, `%rand-zero-modulus`), matching `/lib/fixed`'s style. No NaN-returning "soft" errors for parameter mistakes — parameters are programmer input, not data.
- Float-valued *data* edge cases (e.g., a weight list containing NaN) crash rather than propagate.
- All arms total over their asserted domain: no infinite loops. Rejection loops must have a documented, astronomically-bounded expected iteration count; the Saloon counter-window variant gets an explicit crash on window exhaustion.

## 10. Testing

`/tests/lib/rand.hoon`, runnable with `-test %/tests/lib/rand ~`:

1. **Known-answer tests** (blocking): Philox KAT vectors (Random123), SplitMix64 first-N, PCG64 demo output, Lemire `+below` spot values against a reference trace.
2. **Determinism/replay**: same seed → identical sequence; `+fill-uniform` sequential Hoon equals per-element counter formula.
3. **Unbiasedness**: `+below` over small `n` (e.g., 3, 6, 7) for 10k draws — chi-square statistic below a fixed threshold (hardcode the threshold and the expected counts; this is a smoke test, not TestU01).
4. **Distribution moments**: sample mean/variance of normal, exponential, gamma within `+is-close:rd` tolerances at 50k draws with a FIXED seed (so it's a regression test, not a flaky statistical test — record expected values as constants).
5. **Shuffle**: permutation property (same multiset), and uniformity smoke test over all 6 permutations of a 3-list.
6. **Offline** (README instruction, not Hoon): pipe jet output into PractRand/dieharder. Hoon-side statistical testing beyond smoke tests is out of scope.

## 11. Jetting plan (follow-up milestone, spec'd now)

- Register under the existing `%non` jet chapter pattern.
- Jets for: `+block` (Philox), `+next` (SplitMix, PCG), `+below`, and the Saloon `+fill-*` arms (the ones that matter for throughput; the fill jets use per-element counters and OpenMP/threads freely because the spec fixed the bit-exact per-element mapping).
- Distribution transforms need no dedicated jets at first: they decompose into engine draws + `/lib/math` calls, both already jetted.
- Parity harness: like `libmath/tools/rq_check.c`, a `librand/tools/rand_check.c` that runs the C reference (Random123, PCG reference, Vigna's SplitMix) against dojo output over a seed sweep.

## 12. Non-float output types (`++uni` and `++dist` adapters)

Adapters over the same engines; no new engine work. Ordered by difficulty.

### 12.1 Fixed-point (`/lib/fixed`)

The fixed lattice is uniform, so uniform generation is **exact by construction** — stronger than the float case, no rounding argument required.

- `++fixed  |=([r=rng p=prec] [@ rng])` — full-range uniform: draw `N = (wid p)` raw bits, return as the two's-complement pattern.
- `++fixed-unit  |=([r=rng p=prec] [@ rng])` — uniform `[0,1)`: draw `b.p` bits as the fraction field, integer field zero.
- `++fixed-between` — signed range via the twoc adapter (12.2) at width N, then reinterpret.
- Distributions: sample at `@rd`, quantize round-to-nearest (extend `/lib/fixed` with a `+from-rd` mirroring its `+from-rs` if absent). Single quantization step; document it.

### 12.2 Twoc (`/lib/twoc`)

- `++twoc-full  |=([r=rng w=@] [@ rng])` — w raw bits, interpreted as the width-w two's-complement pattern. Done.
- `++twoc-between  |=([r=rng w=@ a=@ b=@] [@ rng])` — a,b as width-w patterns, `a <= b` in twoc order (crash otherwise). Bias to unsigned (subtract a via twoc `+sub`, reinterpret as unsigned offset), Lemire `+below` on span `b - a + 1`, bias back via twoc `+add`.
- **JET WARNING (put this in a source comment):** the span of a width-64 range needs 65 bits. Hoon `@` doesn't care; a C jet using `uint64_t` will overflow. The jet must special-case span = 2^w or compute the span in 128-bit.
- Since `+twid` is keyed on arbitrary widths, this adapter is the primitive `++fixed-between` delegates to.

### 12.3 Complex (`/lib/complex`)

Component-wise over the packed representation; one arm set per width door (`cs`/`cd` first, matching that library's ship order).

- `++cuniform` — uniform over the unit square: two `[0,1)` component draws, `+pak`.
- `++normal-parts` — iid N(0,1) real and imaginary components.
- `++cnormal` — standard circularly-symmetric complex normal CN(0,1): components N(0, 1/2), i.e. `normal-parts` scaled by `invsqt2`. **Both exist and both are named unambiguously** — "complex Gaussian" means different things to signal-processing and statistics users, and conflating them is a guaranteed bug report.
- `++on-circle` — uniform on the unit circle: `theta = tau * u`, components via `/lib/math` `+cos`/`+sin` (bit-specified kernels, jet parity holds).
- `++in-disk` — uniform on the unit disk: rejection from the square (variable-consumption arm, document per section 5 policy).

### 12.4 Unum / posits (`/lib/unum`)

Two inequivalent uniform semantics. Ship both, named so they cannot be confused:

- `++posit-lattice  |=(r=rng [@ rng])` per width door — uniform over **bit patterns** excluding NaR: draw n bits, redraw on the NaR pattern (single-pattern rejection, expected 1 + 2^-n draws). The induced *value* distribution is approximately log-uniform (posits taper: dense near ±1, sparse at extremes). This is a fuzzing/property-testing primitive, not a statistics primitive; the docstring says so.
- `++posit-unit  |=(r=rng [@ rng])` — uniform **value** on `[0,1)`, exactly: draw k raw bits `u`, form the g-layer dyadic `[%p %.y (dif:si --0 (sun:si k)) u]` (`%z` when u = 0), encode via `+bit` — one RNE round. Correctness condition: every posit rounding boundary in `[0,1)` at width n is dyadic with exponent >= -(4(n-2)+1) (minpos = 2^-4(n-2)), so any k above that bound makes the output *exactly* the round-to-nearest image of uniform measure — not approximately. Fix **k = 32 (posit8), 64 (posit16), 128 (posit32)**. No float intermediary, no double rounding; this construction exists only because `@` is arbitrary-precision and `+bit` is a single rounding step.
  - Documented standard-conformant wrinkle: posits never round nonzero to zero (`+bit` saturates underflow to minpos, per the 2022 standard), so P(zero) = 2^-k exactly and the mass below ~minpos lands on minpos. Not a bug; write the comment so nobody "fixes" it.
- Distributions: sample at `@rd`, convert via the existing `+from-rd`. Double rounding is harmless at p8/16/32 (float64's 52 fraction bits strictly dominate posit32's max 27-bit fraction) — conveniently the same widths where `/lib/unum`'s transcendentals are verified. **p64/p128 inherit that file's existing caveat verbatim: treat as unverified until the oracle sweep is extended.**
- Tests: verify `++posit-unit` at posit8 against an mpmath oracle that computes each posit8 value's exact rounding-interval measure; chi-square the empirical counts at a fixed seed. Lives alongside `unum_cheb_check.py` in `librand/tools/`.
- NEXT-STEPS entry: quire-accumulated Monte Carlo — `+fdp` makes sample *sums* singly-rounded, a capability hardware floats don't have. Design note only; out of scope for v1.

## 13. Milestones for Sonnet

1. `++split-mix` + `++seed` + KATs. (Small; validates the test harness.)
2. `++philox` + KATs. `++fork` across all engines + path-sensitivity tests. `++gen` facade (thin; lands here so later milestones can use it in Examples blocks).
3. `++uni` complete (bits/below/between/floats) + bias tests.
4. `++pcg` + jump.
5. `++dist` at `@rd`: normal, expon, gamma, beta, bernoulli, geometric. Moment tests.
6. `++sample`: shuffle/permutation/choice/alias. `++dist` categorical/poisson/binomial/dirichlet, `@rs` variants.
7. Non-float adapters (section 12): twoc + fixed (one sitting — twoc is the primitive), then complex, then posit lattice/unit + mpmath oracle test.
8. Saloon `+rand-ray` with the counter-window design + replay-equality tests.
9. README + NEXT-STEPS (ziggurat, BTPE, buffered Philox, posit rays, quire Monte Carlo note, `@rh`/`@rq` distributions).

Each milestone lands with tests green before the next starts. If any reference constant in this spec conflicts with the cited reference implementation, **the reference implementation wins** — flag the discrepancy rather than silently choosing.
