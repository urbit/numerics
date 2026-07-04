# Measured max-ULP results: `/lib/math` and `/lib/unum` Chebyshev transcendentals

Generated 2026-07-03 by running `cheb_check.py` and `unum_cheb_check.py`'s full
function/precision grid (extended this session to cover `@rh`, full `@rq`,
`sqt`/`cbrt`/`pow`/`pow-n`/`log-2`/`log-10`/`atan2` for math, and
`cbrt`/`pow`/`pow-n`/`tan`/`log-2`/`log-10` for unum). Raw run output:
[`ULP-RESULTS-raw.txt`](ULP-RESULTS-raw.txt) — the full `@rd`/@rs/@rh/@rq
math grid, the unum `pow`/`pow-n`/`cbrt` additions, and (appended) a re-run
of unum's `log2log10`/`trig`/`exp`/`log`/`atan`/`ainv` for the record (`sqt`
was verified via a small standalone script against `posit_check.py`'s
`my_sqrt`, not one of the file's own `check_*` keys — see the snippet at the
end of the raw log).

Methodology per the paper text: `@rh` is exhaustive over the full 2^16 space;
`@rs`/`@rd`/`@rq` are sampled (broad range + edge-weighted near breakpoints/
powers-of-2), since exhaustive testing is infeasible at those widths. `mpmath`
at high precision is NOT the oracle — each `check_*` function's own exact
`Fraction`/big-int arithmetic is (mirrors `cheb_check.py`'s own docstring:
mpmath independently verifies the algorithm-of-record, it doesn't define it).

## Table A: `/lib/math` (IEEE 754), max observed ULP

| Function | @rh | @rs | @rd | @rq |
|---|---|---|---|---|
| exp | 1.002 (exh.) | 0.860 | 0.994 | 0.806 |
| log | 0.841 (exh.) | 0.613 | 0.986 | 0.617 |
| log-2 | 1.566 (exh.) | 0.901 | 1.210 | 0.600 |
| log-10 | 1.370 (exh.) | 1.303 | 1.459 | 1.103 |
| sin | 0.505 (exh.)† | 0.646 | 0.727 | 0.726 |
| cos | 0.531 (exh.)† | 0.534 | 0.760 | 0.718 |
| tan | 1.637 (exh.)†,∥ | 0.939∥ | 0.757‡ | 1.668∥ |
| atan | 0.736 (exh.) | 0.575 | 0.713 | 0.533 |
| atan2 | 1.800 (sampled) | 1.159 | 1.362 | 1.350 |
| asin | 1.010 (exh.) | 0.849 | 0.837 | 0.765 |
| acos | 0.996 (exh.) | 0.980 | 0.887 | 0.856 |
| sqt | 0.500 (exh.) | 0.500 | 0.500 | 0.500 |
| cbrt | 6.700 (exh.)§ | 2.623§ | 4.394§ | 4.288§ |
| pow | 31.0 (sampled)¶ | 48.0 | 133.7 | 60.7 |

† **FIXED (2026-07-03), numbers below are post-fix.** `@rh` sin/cos/tan
previously did quarter-turn reduction entirely in native `@rh` arithmetic —
every constant piece and intermediate product was itself only an 11-bit-
mantissa value. Half precision only guarantees exact integers up to 2048;
once the reduction quotient `q=round(|x|*2/pi)` exceeded that (`|x|` past
~500), `q`'s own rounding error (up to ~1) got multiplied by `pi/2`'s hi
part, injecting a multi-radian error into the reduced remainder. The
**pre-fix** numbers were catastrophic: max **437,291 ULP (sin)**,
**1,818,711 ULP (cos)**, **21.6 billion ULP (tan)** — not previously
documented, found while writing the paper's accuracy section. **Fixed** by
widening to `@rs` (24-bit mantissa, exact `q` up to 2^24 — vastly beyond
`@rh`'s entire dynamic range) via the existing-but-previously-unused
`+widen-hs`/`+narrow-sh` primitives — exactly the architecture `+widen-hs`'s
own docstring already claimed but was never wired up for trig. Verified
exhaustively in the Python oracle AND live on a fresh ship (bit-for-bit
match at the old blowup points; full `math-rh`/`math-trig`/`math-tan`/
`math-derived`/`math-atan`/`math-ainv` test suites green). The C jet was
also fixed to match (same widen/narrow pattern via SoftFloat's own
`f16_to_f32`/`f32_to_f16`) and re-verified live with jets enabled — every
case now runs in 50-65µs (jetted), not the ~8ms interpreted timing, and
still bit-exact with the Hoon spec. See `libmath/NEXT-STEPS.md` for the
full writeup; the Hoon/jet mismatch this fix could have introduced is
closed, not open.

‡ `@rd` tan has a dedicated fdlibm kernel (`__kernel_tan`, 0.757 ULP) *and* a
sin/cos-ratio fallback path (1.934 ULP) — the table reports the dedicated
kernel's number, which is what `+tan` actually calls.

§ `cbrt` at all precisions is `sign(x)*exp(log|x|/3)`, a composed algorithm,
not independently minimax-fit — its error is inherited from `exp`/`log`
composition, not directly controlled. **Not faithfully rounded** (2.6–6.7
ULP observed) despite the paper text's blanket "faithfully rounded" framing;
this section's language should be corrected (see below).

¶ `pow` (general `x^y`, not the integer-exponent `pow-n` fast path) is
`exp(y*log x)`, likewise composed and **not faithfully rounded** — 31–134
ULP observed depending on precision, worst near extreme `(x,n)` pairs (e.g.
`x=0.01, n=-20` at `@rd`). `pow-n`'s dedicated integer-exponent path is much
tighter (spot-checked at `@rd`: 0.5 ULP at n=2, 1.2–1.9 ULP at n=3/5, 11.5 ULP
at n=-3 — the negative-exponent case routes through division and inherits
more error).

∥ **A real bug in the oracle itself, caught by a reviewer — and the trail led
to a genuine `@rs` accuracy improvement.** Originally, `tan` at `@rh`,
`@rs`, and `@rq` was all `(div (sin x) (cos x))` in `math.hoon` — composed,
no dedicated kernel — only `@rd` had one (`+rd-tan`, a genuine ported
fdlibm `__kernel_tan`). `cheb_check.py`'s `check_tan_rs()` was nonetheless
testing a *hypothetical* dedicated f32 tangent kernel (`tan_f32`/`ktan32`,
mirroring `@rd`'s real one almost coefficient-for-coefficient) that had
apparently been drafted at some point but never actually shipped in Hoon —
so the originally-reported 9.675 ULP was measuring an algorithm that didn't
exist in the codebase, not `@rs`'s actual `sin`/`cos` ratio (1.216 ULP, in
line with its `@rh`/`@rd`/`@rq` neighbors at the time).

Investigating *why* that abandoned draft scored worse than the ratio it was
meant to replace (a dedicated kernel should never lose to a naive
composition) turned up the real story: the draft multiplied its dominant
linear term through the whole polynomial product instead of adding it
**last**, the way fdlibm's own `kernel_tan` structures the computation.
Fixing just that ordering — reusing the draft's Chebyshev-fit-quality
coefficients — brought `@rs`'s *actual, native* dedicated kernel down to
**0.939 ULP**, beating the composed ratio outright. This is now shipped:
`+rs-tan` in `math.hoon`, the C jet `_rs_tan`, verified bit-exact on a live
ship with jets enabled. `@rh` and `@rq` were also investigated (widening
the same approach) but a dedicated kernel measurably regresses `@rh`
(~2.25 ULP native vs. 1.64 composed — 11-bit mantissa, not enough headroom)
and `@rq` needs a fundamentally different exact-arithmetic technique (its
series coefficients grow rather than shrink at high degree, destabilizing
native chained rounding) — both are deliberately left composed. Full
writeup in `libmath/NEXT-STEPS.md`.

Separately, audited every other `check_*_rs`/`_rq`/`_rh` function for the
same failure mode (testing an aspirational algorithm instead of the shipped
one) — `pow`, `cbrt`, `atan2`, and `sqt` at all four precisions, plus the
(unaffected) `tan` at `@rd`/`@rh`/`@rq`, all correctly self-describe as
composed and their Python bodies genuinely call the same composition
math.hoon does. The `tan_rs` oracle bug was an isolated incident, not a
systemic problem.

## Table B: `/lib/unum` (2022-standard posits), max observed ULP

| Function | @rpb (posit8) | @rph (posit16) | @rps (posit32) |
|---|---|---|---|
| exp | 0 (1,143 tested) | 0 | 0 |
| log | 0 (13,984 tested) | 0 | 0 |
| log-2 | 0 (13,984 tested) | 0 | 0 |
| log-10 | 0 (13,984 tested) | 0 | 0 |
| sin | 0 (21,980 tested) | 0 | 0 |
| cos | 0 (21,980 tested) | 0 | 0 |
| tan | 0 (21,980 tested)† | 0 | 0 |
| atan | 0 (11,970 tested) | 0 | 0 |
| asin | 1 (286 tested) | 1 | 3 (1 non-faithful) |
| acos | 1 (286 tested) | **7** (6 non-faithful) | 4 (4 non-faithful) |
| cbrt | 1 (572 tested, x>0 only)‡ | 2 (15 non-faithful) | 2 (8 non-faithful) |
| pow | 3 (1,682 tested, 25 non-faithful) | 4 (40 non-faithful) | 4 (46 non-faithful) |
| pow-n (worst over n=2,3,5,7) | 2 (n=7, 16/600 non-faithful) | 2 (n=5/7) | 2 (n=7, 2/600 non-faithful) |

† `tan` had no dedicated oracle check before this session (only `sin`/`cos`
were covered); added by directly exercising the module's existing `tan_g`
helper (already present, just never wired into a `check_*` function) over
the same grid as `check_trig()`. Result: **correctly rounded (0 ULP)**,
better than the README's previous "faithful (a few ULP)" characterization —
`+tan`'s single-final-rounding `gdiv` on exact `sc`-derived intermediates
turns out to lose no accuracy at these widths, so that characterization
should be tightened.

‡ `/lib/unum`'s `cbrt` returns NaR for `x<0` (unlike `/lib/math`'s
sign-extraction trick), so the sweep is restricted to `x>0`.

`log-2`/`log-10` also had **no dedicated ULP sweep at all** before this
session — `check_log2log10()` printed a handful of bit patterns and stopped;
it now runs the same broad + near-power-of-2 grid as `check_log()` (13,984
points). Result: correctly rounded (0 ULP) at all three widths, as expected
given it reuses `+lr`'s already-correctly-rounded mantissa/exponent split.

## Process notes (bugs caught while producing this report)

1. **A real bug in this test harness, caught and fixed**: the new `@rh`
   section in `cheb_check.py` originally defined a module-level constant
   `PI_H` (half-precision π, for `atan2_f16`) that silently **shadowed** the
   pre-existing `PI_H` (`@rd`'s full-precision π, used by the *original,
   unmodified* `acos_f64`) — Python has one flat module namespace, so the
   later definition wins for every subsequent call, not just calls textually
   below it. This corrupted the long-standing `@rd` `acos` check to report a
   fabricated **2.18 trillion ULP** "error" (`(true acos - computed)` blew up
   because `acos_f64` was silently using `PI_H ≈ 3.140625` instead of
   `3.14159265358979...`). Caught by cross-checking against `math.hoon`'s own
   doccord examples, not by inspection — a reminder that even a "just
   printing extra numbers" script edit can silently corrupt unrelated,
   already-working checks via global name collision. Fixed by renaming to
   `PI_RH`; a full-file scan afterwards confirmed no other top-level name
   collisions between new and pre-existing constants.
2. **An accidental O(2³²) sweep, caught by an 18-minute hang.** The first
   draft of `check_pow_rh`/`check_atan2_rh` (both two-argument functions)
   swept two independent ranges of `2^16` half-precision bit patterns each
   with a small stride — arithmetically still `~65536² ⁄ stride²`  pairs, e.g.
   `atan2`'s original `range(0,65536,3) × range(0,65536,37)` is **~38.7
   million** `(y,x)` pairs, each doing an `mpmath` call. This alone accounted
   for essentially all of an 18+ CPU-minute run that should have taken ~3
   minutes total. Fixed by bounding both to explicit ~400-point grids
   (consistent with the two-argument `@rd`/`@rs`/`@rq` sweeps elsewhere,
   which already use bounded grids, not raw bit-pattern strides). Flagging
   in case this pattern gets copied elsewhere in future two-argument `@rh`
   checks.
3. **Several existing `math.hoon` doccord example values appear stale**,
   found while cross-checking this session's simulation against them (used
   as a sanity check that the Python composition matches the real Hoon
   algorithm): the default-tolerance examples for `@rd` `pow(2,3.5)`,
   `cbrt(2)`, and `log-10(0.1)` don't match either the mathematically true
   value or their own `[%z 1e-15]` tight-tolerance sibling example a few
   lines below in the same docstring — while this session's independently
   bit-validated simulation (`exp`/`log` confirmed bit-identical to
   `math.hoon`'s baked hex before doing anything else) matches the
   tight-tolerance examples almost exactly. `@rq`'s `cbrt(2)` doccord
   (`1.2598919398737178526805575821133312`) is off from the true `2^(1/3)`
   (`1.2599210498948731647672...`) by ~3×10⁻⁵ relative — both its default
   and `[%z 1e-10]` variants show the same wrong-ish value, suggesting a
   copy-paste error in the docstring rather than independent staleness. Not
   fixed here (out of scope for an accuracy-*measurement* pass) but flagged
   for a follow-up doc pass.

## Headline comparison

Every unum transcendental that has a **dedicated, non-composed** kernel
(exp/log/log-2/log-10/sin/cos/tan/atan) is **exactly correctly rounded (0
ULP)** at all three tested widths — a stronger result than `/lib/math`
achieves at any precision (math's best case is ~0.5–1 ULP faithful, not
correctly-rounded 0 ULP), because unum's arbitrary-precision Hoon atoms let
the reduction/kernel arithmetic run at whatever fixed-point width the
worst case needs, with no hi/lo-split compromise forced by a fixed hardware
register width. The **composed** functions in both libraries (`pow`, `cbrt`,
and unum's `pow-n` for exponents beyond the exact/cheap cases) lose this
property in both libraries alike, landing in the 1–130 ULP range depending on
precision and how many chained roundings the composition involves — this is
an inherent cost of composition, not a library-specific weakness.
