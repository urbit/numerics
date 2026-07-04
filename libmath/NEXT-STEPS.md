# `/lib/unum` — next steps

Status as of 2026-06-01.  The posit scalar layer (PR #13, merged) plus the
follow-up domain fixes (PR #38) implement a complete, SoftPosit-verified
`posit<n,2>` library: encode/decode, arithmetic, sqrt/round/fma, integer and
IEEE-754 conversion, the quire + fused dot product, and naive transcendentals,
at widths posit8/16/32/64/128.  This file tracks what is deliberately *not* yet
done, roughly in dependency order.

## Roadmap — 2026-07-03 (Chebyshev transcendentals)

Item 4 below is **done**: `/lib/unum`'s naive Taylor/AGM transcendentals are
replaced with range-reduced Chebyshev-minimax/exact-Taylor kernels, mirroring
`/lib/math`'s own Chebyshev rewrite but exploiting Hoon's arbitrary-precision
`@` atoms (no hi/lo constant splitting needed — decode once via `+sea`, exact
bignum arithmetic via new shared `+gmul`/`+gadd`/`+gneg`/`+gsub`/`+gdiv`/
`+gpoly`/`+g-round`/`+glt` helpers, round once via `+bit`). `exp`/`log`/
`log-2`/`log-10`/`sin`/`cos`/`atan` are now correctly rounded (0 ULP vs
mpmath) at posit8/16/32; `tan`/`asin`/`acos` are faithful (composed from the
new `atan` plus existing correctly-rounded ops, not a dedicated rational
kernel — worst case observed is 7 ULP for `acos` at posit16, not merely "a
few"). Also fixed a pre-existing `+acos` bug (wrong quadrant for x<0).

**Scope gap, not yet closed**: these arms live in the generic `bloq`-
parameterized `+pp` core with no width guard, so `+rpd`/`+rpq` (posit64/128)
now run this same new code too — but accuracy there is UNVERIFIED (only
posit8/16/32 were checked against the oracle; `WBITS=128` was sized for
posit32's worst case, with no margin proven sufficient at posit128). This
supersedes item 3's framing below ("rpd/rpq jets return u3_none, fall back to
Hoon" is still true for the JET, but the Hoon itself is no longer the old
naive series either — it's the new, unverified-at-that-width Chebyshev code).

SoftUnum (the C jet twin) is ported to match, using **GMP** (`mpz_t`) for the
transcendental kernels' exact g-layer arithmetic — the existing fixed-512-bit
`wide_t` can't replicate Hoon's never-truncate-until-final-round approach
without its own precision-margin proof, and GMP already has precedent here
(`/lib/twoc`'s jet). Along the way, found and fixed a real bug: GMP's
allocator is hooked to vere's `u3a` per-road heap (`mp_set_memory_functions`
in `manage.c`'s `u3m_init()`), so any C state cached in a `static` variable
across separate jet calls silently corrupts once its allocating road's heap
watermark resets on the next top-level call — no jet in vere caches raw heap
state across calls (checked `/lib/twoc` too), so SoftUnum's transcendental
kernels now recompute their shared constants fresh every call and clear them
before returning, same convention as everywhere else.

Benchmarked at `benchmark/results/2026-07-03/unum-chebyshev/`: interpreted
transcendentals are **33–77× faster** (avoiding the naive series' repeated
posit decode/encode per term); jetted is roughly a wash except `atan`, now
**25–46× faster jetted at posit16/32** — the old AGM-based atan hit the
`wide_t` bit-by-bit-`isqt` bottleneck (item 2 below) on every iteration, the
new fdlibm-breakpoint atan has no `sqrt` in it at all, sidestepping that
bottleneck as a side effect of an accuracy-motivated rewrite.

PRs: numerics #71 (Hoon rewrite), SoftUnum `930fe6d` (GMP-backed C port,
pushed to master), vere #1046 (updated in place — re-pinned `ext/softunum`,
linked GMP via the already-vendored `ext/gmp`, `jets/i/unum.c` needed no
changes).

Still open from this item: quire-accumulated summation inside the polynomial
evaluation (the more ambitious original idea here) wasn't needed — exact
bignum Hoon-atom arithmetic already gives single-rounding correctness
without it. Not revisited.

## `/lib/math` — `@rh` sin/cos/tan large-argument bug, FIXED (2026-07-03)

`math.hoon`'s `@rh` (half-precision) `+rh-trig` engine did quarter-turn range
reduction (`q=round(|x|*2/pi)`, `r=|x|-q*(pi/2)`, 3-part `pi/2` split) entirely
in native `@rh` arithmetic — every constant piece and every intermediate
product was itself only an 11-bit-mantissa half-precision value. Half
precision only guarantees exact integer representation up to 2048; once
`q` exceeds that (`|x|` past ~500), `qf=(sun q)` silently rounds to the
nearest even `@rh`-representable integer, and that error (up to ~1) gets
multiplied by `pi/2`'s hi part (1.5), injecting a multi-radian error into the
reduced remainder. Measured impact: max ULP error grew from faithful (<1
ULP) below `|x|~500` to **437,291 ULP (sin)**, **1,818,711 ULP (cos)**, and
**21.6 billion ULP (tan)** by `|x|~3000-43000` — well within `@rh`'s normal
range (max ~65504), not an edge case. Not previously documented; found while
writing the paper's accuracy section (`doc/` USTJ manuscript).

**Fixed** by widening to `@rs` (single precision, 24-bit mantissa — exact
integer `q` up to 2^24, vastly beyond anything `@rh`'s dynamic range can
demand) via the existing-but-previously-unused `+widen-hs` (exact) and
narrowing back via `+narrow-sh` (correctly-rounded RNE) — this is exactly the
architecture `+widen-hs`'s own docstring already claimed ("The @rh
transcendentals compute in the (more precise) @rs door"), just not actually
wired up for trig before now. The old native-`@rh` `+rh-trig` engine
(`sc`/`cc`/`neg`/`ksin`/`kcos`/`trig-fin`) is now dead code, removed.
Verified: Python oracle (`cheb_check.py`) exhaustively over the full `@rh`
domain gives max 0.505 ULP (sin), 0.531 ULP (cos), 1.637 ULP (tan) — matches
live on-ship behavior exactly (spot-checked at the old blowup points, plus
the full `math-rh`/`math-trig`/`math-tan`/`math-derived`/`math-atan`/
`math-ainv` test suites, all green). New regression tests added to
`tests/lib/math-rh.hoon` covering the old blowup points.

**Hoon/jet mismatch: CLOSED (2026-07-03).** `libmath/vere/noun/jets/i/math.c`
and `vere64`'s copy previously still implemented the OLD, superseded
native-`@rh` algorithm line-for-line (jets for `sin`/`cos`/`tan`\@`@rh` are
registered `no_hashes`, i.e. matched by name only), so a jetted host was
silently computing a **different, wrong** answer than the corrected Hoon
specification for large `|x|` — a real jet mismatch by this project's own
architectural discipline (Section 3 of the USTJ paper). Discovered while
ship-verifying the Hoon fix: testing the Hoon change required temporarily
stripping the `~/  %sin`/`~/  %cos`/`~/  %tan` jet hints to observe pure-Hoon
behavior, since the stale jet otherwise silently overrode every test — a
reminder that `no_hashes` jets are invisible landmines for any future
algorithm change to these arms.

**Fixed** by porting the same widen-@rs/narrow-@rh pattern to C: `_rh_sin`/
`_rh_cos` now call SoftFloat's own `f16_to_f32` (exact) and `f32_to_f16`
(correctly rounded RNE) around the existing, already-correct `_rs_sin`/
`_rs_cos`, removing the old native-`@rh` kernel
(`_rh_ksin`/`_rh_kcos`/`_rh_trigfin`) entirely — both conversion primitives
were already compiled into the linked SoftFloat library, so no new code
beyond the jet functions themselves was needed. Applied to the real vere
source (`pkg/noun/jets/i/math.c`) and mirrored to both
`libmath/vere/noun/jets/i/math.c` and `libmath/vere64/noun/jets/i/math.c`.
Rebuilt vere (`zig build`, aarch64-macos-none target) and verified on a
freshly booted ship: the full `tests/lib/math-rh.hoon` suite (including the
large-`|x|` regression cases added for this fix) passes with jets **enabled**,
and every case runs in 50-65 microseconds — confirming the new jet fires
(not falling back to interpretation) and is bit-exact with the corrected
Hoon specification, closing the mismatch.

## `/lib/math` — `@rs` dedicated `tan` kernel, DONE; `@rh`/`@rq` deliberately left composed (2026-07-03)

Prompted by a reviewer noticing `@rs` `+tan`'s 9.68 ULP figure (Table
tab:ulp-math) was an unexplained outlier next to its 0.5–1.7 ULP neighbors:
investigation found the 9.68 ULP number was itself an oracle bug (see the
`tan_rs`-mismatch entry above), and that fixing the oracle revealed `@rs`
was on the plain `(div (sin x) (cos x))` composed path (~1.22 ULP) the
whole time — only `@rd` has ever had a genuine dedicated kernel
(`+rd-tan`, a real ported fdlibm `__kernel_tan`). This prompted the
question: why not give every precision a dedicated kernel, for uniformity?

**Investigated all three (`@rh`, `@rs`, `@rq`) via the Python oracle before
touching any Hoon, mirroring the project's established discipline.**
Findings, in order of increasing surprise:

- **`@rs`: a clean, real win — DONE.** A properly Chebyshev-fit low-degree
  (7-coefficient) polynomial for `Q(z)=(tan(r)/r-1)/z`, evaluated in
  genuinely native `@rs` arithmetic (no borrowing a wider precision),
  reaches **0.939 ULP** — beating the composed ratio's 1.22 and even
  `@rd`'s own dedicated kernel (0.757 is close; 0.939 is competitive).  The
  key insight, found by trial: the kernel's dominant linear term (`rhi`)
  must be **added last** (`w2 = rhi + r`, mirroring fdlibm's own structure)
  rather than multiplied through the whole polynomial product — the
  abandoned `@rs` draft this project already had (see the oracle-bug entry
  above, `tan_f32`/`ktan32`, 9.68 ULP measured) got exactly this wrong,
  which is *why* it scored worse than the ratio it was meant to replace,
  not because a native `@rs` kernel is inherently a bad idea. A raw
  16-term exact-Taylor attempt (no minimax-quality fit) also
  underperformed at ~1.5 ULP — degree matters as much as structure once
  enough terms accumulate chained-rounding error. Shipped: `+rs-tan` door
  in `math.hoon` (`+redq`/`+ktan`/`+main`, mirroring `+rd-tan`'s shape),
  the C jet `_rs_tan` (ported to the same algorithm using SoftFloat ops
  directly), mirrored to both `libmath/vere` and `libmath/vere64`. Verified
  on a fresh ship with jets enabled: `tests/lib/math-tan.hoon`'s new
  `test-tan-rs-*` cases plus the full `math-trig`/`math-derived`/
  `math-atan`/`math-ainv`/`math-rh` regression suite all green, every case
  running in ~50µs (jetted).  One test-vector correction along the way: the
  `x=pi/4` exact-tie case (`ax*2/pi` lands precisely on `.5`) exposed a
  genuine rounding discrepancy between `cheb_check.py`'s own
  `reduce_pio2_32` (gave `q=1` at this tie) and Hoon's actual
  round-to-even `+redq` (gives `q=0`) — checked against `mpmath` directly,
  Hoon's `q=0` answer (`0x3f800000`) is the mathematically closer one, so
  the test vector was wrong, not the shipped code; harmless for the
  measured ~0.94 ULP figure since a continuous 200k-point sweep essentially
  never lands exactly on a tie.

- **`@rh`: a dedicated kernel would be a regression — NOT done, by
  design.** Even after fixing the reduction (reusing the same widen-based
  approach that fixed `@rh`'s sin/cos), the best native-`@rh` kernel found
  plateaus around **2.25 ULP** — *worse* than the composed ratio's 1.64.
  `@rh`'s 11-bit mantissa doesn't leave enough headroom for a multi-step
  polynomial-plus-reciprocal kernel to beat two already-accurate `+sin`/
  `+cos` calls and one division. Left composed.

- **`@rq`: a dedicated kernel is possible but needs a fundamentally
  different technique — NOT done, deferred.** `tan`'s series coefficients
  *grow* in magnitude at high degree (unlike `sin`/`cos`'s shrinking
  factorial-decay coefficients), and `@rq` needs 45-60+ terms to converge
  at 112-bit precision. Native chained `@rq` arithmetic over that many
  growing-magnitude terms is catastrophically unstable (errors observed in
  the billions of ULP, worsening — not improving — with more terms, a
  classic symptom of the wrong technique rather than insufficient degree).
  The only approach that worked (0.498 ULP) computed the kernel using
  exact/unrounded arithmetic internally, rounding to `@rq` once at the very
  end — mirroring `/lib/unum`'s own g-layer philosophy, but a genuinely
  different internal technique than `@rd`'s or the new `@rs`'s native
  chained-arithmetic style, and not attempted in Hoon here. Left composed
  (~1.67 ULP, already faithful) as a known, explicitly-scoped gap rather
  than a silently-incomplete uniformity story.

**Net effect**: "the same algorithm everywhere" turned out not to be the
right frame — `@rd`'s specific hi/lo-split minimax technique doesn't scale
cleanly to either a much coarser (`@rh`) or much wider (`@rq`) precision.
The honest per-precision picture is now: `@rs` and `@rd` have dedicated
kernels (0.94 and 0.76 ULP), `@rh` and `@rq` use the composed ratio (1.64
and 1.67 ULP) because that's what each precision's own constraints
actually support best today.

## Roadmap — 2026-06-28 (post-jets)

The jet effort (old §5) is **done**: **SoftUnum** (`sigilante/SoftUnum`, the
pure-integer C twin of `/lib/unum`, verified four ways) is vendored into vere
(`ext/softunum`) and the full `/lib/unum` surface (54 arms, posit8/16/32) is
jetted under `non/unum`, fired bit-exact on a hoon-135 fakezod, and benchmarked
(arithmetic ~180–295×, transcendentals ~760–2000×, `fdp` ~760–2300× over
interpreted).  PRs: numerics #65 (hints + reference + valids appendix +
benchmarks); vere #1046 (draft, stacked on the math-jets PR #1044).

Prioritized next steps:

0.  **Land the PRs.**  numerics #65 → merge.  vere #1046 is gated on #1044
    (math jets) and #1022 (lagoon SoftBLAS) merging; then re-target #1046 to
    `develop` and mark ready.  Registration is hoon-135 only (lowest kelvin).

1.  **Posit linear algebra — the Saloon/Lagoon payoff (highest value).**
    `%unum` is already wired into Lagoon (PR #42) and Saloon eig (PR #58), so
    Saloon-over-posits works and now inherits the *scalar* jet speedup (its
    inner `mmul`/`dot`/`sqt`/`add` hit the jetted unum ops).  The remaining win
    is **array-level**: a SoftUnum-backed **posit GEMM/dot using the quire**
    (the analog of SoftBLAS), vendored like SoftBLAS, with a Lagoon `%unum`
    `dot`/`mmul` jet dispatching to it -- one exact-accumulated C call instead
    of a Hoon loop over per-element `fdp`.  This is the headline posit feature:
    exact dot product -> matmul with no error accumulation.  Saloon
    decompositions inherit it for free.

2.  **Perf gaps the benchmark surfaced.**
    - posit16/32 `sqt` is still slow even jetted (jetted `sqt:rps` ~17 us, was
      ~17 us before too — unchanged code): the 512-bit `wide_t` **bit-by-bit
      `isqt`** is the bottleneck.  Replace with a faster wide sqrt (Newton +
      wide divmod, or an `__int128` seed then refine).  `atan`'s OWN instance
      of this bottleneck is gone as of the 2026-07-03 Chebyshev rewrite (the
      new fdlibm-breakpoint atan has no `sqrt` in it) — `sqt` itself, and
      anything else that still calls it (e.g. `asin`/`acos`'s `1-x^2` step),
      still hits the slow path.
    - posit16 over-uses the 512-bit `wide_t` (its quire is 256-bit, arithmetic
      fits `__int128`); a tighter p16 path would shave the common arithmetic.

3.  **Coverage: posit64/128.**  SoftUnum does 8/16/32; `rpd`/`rpq` jets return
    `u3_none` (fall back to Hoon).  The `wide_t` machinery already exists --
    extend SoftUnum to 64/128 (1024/2048-bit quire) to jet the last two widths.
    No external oracle (cerlane `pX2` caps at 32); verify vs the from-scratch
    reference + the Hoon.

4.  **DONE (2026-07-03) — Accuracy: range-reduced transcendentals** (old §2).
    See the "Roadmap — 2026-07-03" section above for the full writeup.
    Quire-accumulation wasn't needed (exact bignum Hoon-atom arithmetic
    already gives single-rounding correctness).

5.  **Standard-name alias layer** (old §1): the 2022-standard public names
    (`addition`/`subtraction`/`sin-pi`/`compound`/`hypot`/`arctan2`/...) over
    the implemented core; mostly renames + thin compositions.

6.  **Valids** (the Type-III interval class) -- see the appendix below.  Lowest
    priority: new surface, no SoftPosit oracle, loosely standardized.

7.  **Parallel Rust SoftUnum** (for NockApp; not urgent).  C and Rust can't
    share code, so "same behavior" is enforced by a language-neutral
    test-vector corpus generated from the oracle that both reproduce.  The C
    repo reserves a top-level `rust/` seam.

8.  **Upstream SoftUnum -> `urbit/SoftUnum`** (like `urbit/SoftBLAS`) once it
    stabilizes, and re-pin the vere `ext/softunum` tarball there.

9.  **Decimal printer / literal syntax** for the posit auras (`@rpb`..`@rpq`):
    emit the §6.3 minimum significant digits.  Optional local runtime patch,
    decoupled from the library.

## Verification & tooling

- **Oracle**: SoftPosit (vendored C in `src/SoftPosit`, plus the `softposit`
  pip package) is the reference for arithmetic and quire ops.  It has **no
  transcendentals** — use `mpmath` → `convertDoubleToPX2` for those.
  `libmath/tools/posit_check.py` is the offline harness (exhaustive posit8,
  sampled posit16/32).
- **On-ship tests** run on a live ship (we used `~hex` via MCP).  There is no
  local urbit binary; a persistent build failure from the MCP tool is a *real*
  compile error — read the dojo (`tmux capture-pane`) for the actual
  `-need`/`-have`, since the MCP surface only returns a generic failure.
- A full per-width Hoon sweep of all 65k posit8 *pairs* times out the test
  runner; keep heavy exhaustion in the Python harness and keep the on-ship
  suite to round-trips, property checks, and curated spot values.

## 1. Standard-name alias layer  (small, additive)

The README's "Posit Standard Compliance" section is the spec: thin arms over
the implemented core.  Most are renames or one-line compositions:

- `negate`/`addition`/`subtraction`/`multiplication`/`division`,
  `compare-{equal,not-equal,greater,greater-equal,less,less-equal}`,
  `sign` ← `sgn`, `nearest-int`/`ceil`/`floor` ← `rnd`/`cel`/`flr`.
- `next`/`prior` — lexicographic successor/predecessor of the bit pattern
  (`+(p)` / `(dec p)` with NaR/extreme handling); genuinely new, but trivial.
- `*-pi` trig (`sin-pi` = `(sin (mul pi x))`), `*-plus-1`/`*-minus-1`
  elementary, `compound`/`root-n`, `hypot`, `fmm`, `arctan2` — compositions.
- inverse and hyperbolic trig (`arcsin`/`arccos`/`arctan`, `sinh`/`cosh`/`tanh`,
  `arcsinh`/…): new naive series, same caveats as the existing transcendentals.

No new infrastructure; can land as one PR.  Decide whether aliases live in
`unum.hoon` itself or a thin `unum-std.hoon` wrapper.

## 2. Transcendental accuracy  (DONE 2026-07-03)

Was: naive fixed-term Taylor series, accurate only near the expansion point.
Now: range-reduced Chebyshev-minimax/exact-Taylor kernels, correctly rounded
(0 ULP vs mpmath) for `exp`/`log`/`log-2`/`log-10`/`sin`/`cos`/`atan` at
posit8/16/32; faithful for `tan`/`asin`/`acos`. See the "Roadmap — 2026-07-03"
section at the top of this file, `libmath/tools/unum_cheb_check.py` (the
mpmath-verified algorithm-of-record), and numerics PR #71. Quire-accumulated
sums (the original idea below) turned out not to be needed — exact bignum
Hoon-atom arithmetic already gives single-rounding correctness without them.

## 3. Lagoon `%unum` integration  (medium; the high-value item)

This is what makes posits useful for arrays/linear algebra.  Lagoon lives on the
**base desk** (`/lib/lagoon`, `/sur/lagoon`) — there is no separate lagoon desk.

- `sur/lagoon.hoon`: the `+$kind` union already has a commented-out `%unum`
  line — uncomment it.  Reconcile the stale aura note there (`@ruw/@ruh/@rub`)
  with the shipped family `@rpb/@rph/@rps/@rpd` (`bloq` selects width).
- `lib/lagoon.hoon` `fun-scalar`: add a `%unum` branch dispatching per-scalar to
  `add:rpb:unum`, `lth:rpb:unum`, etc., keyed by `bloq` (mirror the existing
  `%i754`/`%uint` branches).  Direct reuse of the implemented arms.
- **Reductions/linalg are the payoff**: route `dot`/`mmul`/`sum`/`cumsum`/
  `trace` over `%unum` arrays through the **quire** (`fdp` / repeated
  `q-mul-add` then one `q-to-p`) so products accumulate exactly and round once —
  exact dot product → matmul with no error accumulation.
- `convert`: posit↔i754 via the any-width `to-r*`/`from-r*` matrix (cross-width
  is fine and intended); posit↔uint/int2 via `sun`/`san`/`toi`.
- Jets come later (see §5).

## 4. Valids  (large; the third unum class)

The interval class (`@rvb/@rvh/@rvs`).  Entirely new surface: an interval is a
pair of posit endpoints with open/closed tags; arithmetic is interval
arithmetic.  Lowest priority; design from Gustafson's Type-III definition.

## 5. Jets  (DONE -- see the 2026-06-28 roadmap above)

**Done.**  Rather than vendor SoftPosit (whose `pX2` caps at 32 bits and has no
transcendentals), we built **SoftUnum** (`sigilante/SoftUnum`) -- the bit-exact
pure-integer C twin of this library -- vendored it into vere (`ext/softunum`),
and jetted the full surface (54 arms, posit8/16/32) under `non/unum`.  Verified
firing bit-exact on a hoon-135 fakezod; benchmarked in `benchmark/results/`.
PRs: numerics #65, vere #1046 (draft on #1044).  Follow-on perf/coverage work
(posit GEMM, wide `isqt`, posit64/128) is in the roadmap above.

## Known minor items

- `fdp` silently truncates to the shorter of its two input lists.  Acceptable
  (zip semantics) but undocumented; add a note or `?>` if a caller needs strict
  equal-length.
- Cleanup the review suggested but we deferred: the eight `to-r*`/`from-r*`
  one-liners could collapse to a single `bloq`-dispatched gate; the per-width
  `test-consts-*` / round-trip tests could be table-driven (cf. the lagoon
  test-by-category convention).  Cosmetic.
- `twoc.hoon` now has corrected `lth`/`lte`/`gte` and `overflow`; the only
  former consumer (`lagoon-old.hoon`) was deleted in PR #38.  A real
  `%int2`-into-lagoon effort would be twoc's first live caller.

---

# Appendix: Valids (Type-III interval unums) — assessment & open questions

Status as of 2026-06-28.  Valids are the third Type-III unum (posits · quires ·
**valids**), the rigorous-interval class.  `/lib/unum` does not implement them;
this appendix records the design assessment so the eventual effort starts from a
shared understanding rather than a blank page.  Lowest priority (no consumer, no
oracle, loosely standardized) — slot it *after* the posit/quire jets.

## What a valid is

A valid is a *guaranteed enclosure* of a real quantity: two posit-like
**endpoints, each tagged with a "ubit" (uncertainty bit)**.  ubit = 0 → the
endpoint is exact (closed); ubit = 1 → "the open interval between this posit and
its neighbour".  So a valid encodes `[lo, hi]` / `(lo, hi)` / half-open — a value
*and* its uncertainty.  Every operation returns the **tightest valid containing
all possible results** (Gustafson's "end of error": a provable bound, not a
single rounded answer).

**The projective twist.**  Valids live on the *projective* real line — the reals
closed into a ring with a single point at the top (the NaR/∞ point, shared with
the posit bit layout).  This lets a valid represent **exterior / wrapping
intervals** ("x < −3 OR x > 5") and unbounded sets ("x > 5") by going the long
way round through infinity.  Consequence: compare / union / intersection are
arithmetic on *arcs of a circle*, not segments of a line — genuinely different
semantics from classical `[lo,hi]` interval arithmetic.

## What it would build on (and the one real new primitive)

Endpoints are posits, so valid arithmetic reuses the SoftUnum / `/lib/unum`
posit core.  The new building block is **directed rounding**: our posit `+bit`
only rounds nearest-even, but valid endpoints must round *outward* — the lower
bound toward −∞, the upper toward +∞ — to keep the enclosure sound.  So we need:

  - a round-toward-±∞ posit encoder (round-down / round-up of an exact value),
  - `+next` / `+prior` (lexicographic successor / predecessor of a posit bit
    pattern — already a TODO in §1 above; `+(p)` / `(dec p)` with NaR/extreme
    handling),
  - the ubit then records residual openness after directed rounding.

Everything else (endpoint add/sub/mul/div) is the existing posit arithmetic.

## Open design questions (pin these BEFORE writing code)

1. **Layout / width.**  Does `@rvb`/`@rvh`/`@rvs` mean a valid *built from two
   `posit{8,16,32}` endpoints*, or Gustafson's "valid⟨2n⟩ = two n-bit
   ubit-posits" packed into a byte/half/single?  The byte/half/single aura names
   mirror posits and are currently ambiguous.  **First decision.**
2. **Projective vs bounded.**  Full Gustafson (wrapping / exterior intervals on
   the projective ring) — the "real" valid — vs a simpler bounded `[lo,hi]`
   interval arithmetic that covers most numeric use with far less machinery.
3. **Set operations & special values.**  intersection / union / complement /
   `is-empty` / `is-everything`, plus the special encodings (empty set,
   all-reals, NaR).
4. **Comparisons.**  Posit ordering == two's-complement of the raw bits; that
   identity does **not** carry to valids — comparison becomes set / containment
   relations, needs its own design.

## Verification posture (harder than posits)

- **No SoftPosit oracle** — SoftPosit is posits + quires only.  Plan mirrors the
  transcendentals: a from-scratch **exact-rational interval reference** (Python
  intervals of `Fraction`s with directed rounding), and optionally **Stillwater
  `Universal`** (C++), which does implement `valid<>` types, as a second oracle.
- **Standardization caveat.**  The 2022 Posit Standard normatively pins posits
  and quires; valids are described in Gustafson's broader work but are *not*
  specified to the same degree.  So this is a clean-room *design* of both the
  Hoon `/lib/unum` arms and the C — NOT a transliteration of a fixed spec, which
  is a different (looser) posture from everything done so far, where `/lib/unum`
  was the bit-exact target.
