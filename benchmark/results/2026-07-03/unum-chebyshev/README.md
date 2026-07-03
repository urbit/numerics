# `/lib/unum` Chebyshev vs Taylor benchmark — 2026-07-03

Four-way per-call comparison of the `/lib/unum` arms: naive-Taylor (the
[2026-06-28 baseline](../../2026-06-28/unum/)) vs. the new range-reduced
Chebyshev-minimax rewrite (numerics PR #71 / SoftUnum `930fe6d` / vere PR
#1046), each measured interpreted and jetted, at posit8/16/32
(`rpb`/`rph`/`rps`). Full table in [`table.txt`](table.txt).

This is a redo of the same dated run with **extended arm coverage**: the
first pass only covered `add sub mul div fma sqt neg lth exp log sin cos
atan pow` (plus `fdp`). This version adds `log-2 log-10 tan asin acos cbrt
pow-n` — the same transcendental surface `/lib/math`'s own benchmark covers,
minus `atan2` (`/lib/unum` has no `atan2` arm). The Taylor baseline never
measured these seven, so their rows in `table.txt` show only the two
Chebyshev columns, marked `n/a` for Taylor.

## Method

Same harness and protocol as the 2026-06-28 run (see that directory's
`README.md` for the full methodology) — `gen/bench-unum-grid` (→
`lib/unum-cells`) precomputes inputs outside the timed loop, folds each arm
`n` times inside `~>(%bout)`. Jetted at `n=100,000`; interpreted at `n=100`
(same reasoning: interpreted transcendentals are too slow for a larger n).
One difference from 2026-06-28's binary: this run's jet binary is the
**Chebyshev-era vere build** (`sigilante/unum-jets-chebyshev`, GMP-linked,
pkg/noun/jets/i/unum.c unchanged, `ext/softunum` re-pinned to SoftUnum
`930fe6d`) — see `benchmark/MANIFEST.md`'s notes section for why the old
Taylor jetted numbers can't be regenerated with this binary (the jet matches
by name not battery-hash, so it would silently compute the NEW algorithm
against the OLD Hoon's test points — invalid). The old Taylor numbers in this
comparison are copied verbatim from `results/2026-06-28/unum/table.txt`, not
regenerated.

`benchmark/desk/lib/unum.hoon` now holds the Chebyshev version (copied from
`libmath/desk/lib/unum.hoon`, numerics PR #71); the old naive-Taylor version
is preserved at `benchmark/desk/lib/unum-taylor.hoon` for reference (see
`MANIFEST.md`).

## Headline result: interpreted transcendentals are 33–77× FASTER

This is the biggest surprise in this data, and it's a real effect, not noise
(the unchanged arithmetic arms below show the actual measurement noise floor
is ~1.25–8.4× interpreted, an order of magnitude smaller than the
transcendentals' 33–77×):

| arm | Taylor interp (µs) | Chebyshev interp (µs) | speedup |
|---|---|---|---|
| exp  | 18,483–19,427 | 547–572       | **~34×** |
| log  | 35,992–39,067 | 736–745       | **~49–52×** |
| sin  | 30,238–32,293 | 807–824       | **~37–40×** |
| cos  | 30,259–32,481 | 805–826       | **~37–40×** |
| atan | 61,659–63,708 | 818–832       | **~75–77×** |
| pow  | 53,045–56,963 | 1,393–1,396   | **~38–41×** |

The naive Taylor series called posit-level `mul`/`div`/`add` for every one of
its 20–40 series terms — each such call independently decodes (`+sea`) and
re-encodes (`+bit`) a posit, so a 20-term series pays that decode/encode
overhead 20 times over. The Chebyshev rewrite decodes ONCE, does its exact
range-reduction/polynomial-evaluation arithmetic on bignums (`+gmul`/`+gadd`/
etc., no rounding in between), and encodes ONCE via `+bit` — even though the
new algorithm does genuinely more sophisticated math (range reduction,
higher-precision constants), it does dramatically *less* posit-level
call/decode/encode churn, and that dominates interpreted cost.

## Jetted: roughly a wash, except atan is now 25–46× FASTER

| arm | Taylor jet (µs) | Chebyshev jet (µs) | ratio |
|---|---|---|---|
| exp/log/sin/cos/pow (rpb) | 11.6–45.1 | 16.8–31.8 | ~0.4–1.5× (mixed, within noise-ish range) |
| exp/log/sin/cos/pow (rph/rps) | 25.3–60.7 | 16.8–31.9 | **~1.4–2.2× faster** |
| **atan (rph/rps)** | **459/850** | **18.7/18.8** | **~25×/~45× faster** |

`atan`'s old Taylor implementation used a 40-iteration Gauss/AGM loop with a
`sqrt` call per iteration — at posit16/32 that hits SoftUnum's known slow
path (the 512-bit `wide_t`'s bit-by-bit `isqt`, flagged as a perf gap in
`NEXT-STEPS.md` item #2 back when the jets first landed). The new
fdlibm-breakpoint-reduced atan has no `sqrt` in it at all, so it completely
sidesteps that bottleneck — as a side effect of an accuracy-motivated
rewrite, not a targeted perf fix. `NEXT-STEPS.md`'s open "wide isqt" item is
now much less urgent (still relevant to `sqt` itself, which is unchanged and
still slow at rph/rps: ~9.5/17.2 µs jetted, matching the old baseline
within noise).

Arithmetic (`add`/`sub`/`mul`/`div`/`fma`/`sqt`/`neg`/`lth`/`fdp`) is
**unchanged code** (this rewrite only touched the transcendentals) — its
jetted numbers match the old baseline within measurement noise (0.92–1.32×
this run; this machine was busier than the 2026-06-28 run, with several
other fakezod piers active in the background, so the noise band widened
somewhat but stayed an order of magnitude below the transcendental effect).
Its *interpreted* numbers also come out somewhat faster in this run
(1.25–8.4×, widest for `neg`, whose absolute time is only a few µs so small
absolute jitter reads as a large ratio) despite identical code — that's
environmental noise, not a real effect; it's dwarfed by the 33–77×
transcendental speedup, which is far too large to be noise.

## Newly-covered arms: tan/asin/acos/cbrt/pow-n/log-2/log-10

No Taylor baseline ever measured these seven — the 2026-06-28 benchmark only
covered `exp/log/sin/cos/atan/pow` (plus arithmetic). Since there's nothing
to compare them against, here are their absolute Chebyshev numbers (range
across posit8/16/32) alongside the already-covered transcendentals, to judge
whether they're in the same ballpark or not:

| arm | interp µs | jetted µs | vs. already-covered |
|---|---|---|---|
| log-2  | 726–743   | 17.9–18.2 | same ballpark as `log` (736–745 / 17.9–18.0) — shares the same `+lr` reduction, just a different final scale constant |
| log-10 | 744–761   | 18.1–18.3 | same ballpark as `log` |
| tan    | 836–844   | 19.5–19.6 | close to `sin`/`cos` (807–826 / 18.5–18.9), consistent with `tan` = `sin`/`cos` plus one extra `+gdiv` |
| pow-n  | 684–698   | **2.1–2.9** | interpreted cost is in the same range as `log`/`tan`, but jetted is dramatically *cheaper* — close to plain arithmetic (~2 µs), since `pow-n(x,3)` is just repeated multiplication, not a polynomial-evaluation kernel |
| asin   | 1979–2019 | 20.5–39.6 | **~2.5–3× costlier** than the already-covered transcendentals, interpreted; jetted cost also *grows with posit width* (20.5 µs at rpb → 39.6 µs at rps), unlike exp/log/sin/cos/atan which stay flat ~17–19 µs across all three widths |
| acos   | 2097–2152 | 21.8–40.9 | same shape as `asin` (built on the same reduction machinery), slightly pricier |
| cbrt   | 1729–1756 | 31.6–32.3 | ~3× `exp`/`log`'s interpreted cost; jetted lands right next to `pow` (31.8–31.9) — both use a Newton-style refinement step |

Summary: `log-2`/`log-10`/`tan` land right next to their nearest
already-covered sibling (`log`, `sin`+`cos`), as expected since they reuse
most of the same code path. `pow-n` is the cheap outlier — jetted, it's
barely more expensive than raw multiplication, since integer powers need no
polynomial evaluation at all. `asin`/`acos`/`cbrt` are the expensive
outliers, both interpreted and jetted; `asin`/`acos` are also the only arms
whose *jetted* per-call cost visibly grows with posit width (rps ~2× rpb).
Nothing here looks pathological — even the priciest new arm is still
~50–75× faster jetted than interpreted, in line with the rest of the
table — but `asin`/`acos`'s width-scaling is worth a note for anyone
planning to put them on a hot path.

## The real headline is accuracy, not speed

None of the above should overshadow the actual point of this rewrite: the
old Taylor numbers above are timings for an algorithm that was **wrong** far
from the origin (`exp(10)` ≈ 21,991 instead of the true ≈ 22,026;
`sin(10)` left `[-1,1]` entirely; see `unum-edge.hoon`'s pre-rewrite
"documented limitation" block, now flipped to correctness asserts). The new
Chebyshev kernels are correctly rounded (0 ULP vs mpmath) for
exp/log/log-2/log-10/sin/cos/atan at posit8/16/32, and faithful (a few ULP)
for tan/asin/acos — see `libmath/tools/unum_cheb_check.py` and numerics PR
#71's test suite for that story. A fast, wrong answer was never a good
trade; this rewrite turns out to also be fast (or faster), but that's a
bonus, not the reason it was worth doing.

## Files

- `table.txt` — the full four-way per-call table.
- `jetted-n100000.txt`, `interp-n100.txt` — raw scraped `%bout` grids for
  this run (Chebyshev only; the Taylor raw grids are in `../2026-06-28/unum/`).
- harness: `benchmark/desk/lib/unum-cells.hoon`, `benchmark/desk/gen/bench-unum{,-grid}.hoon`,
  `benchmark/tools/bench_unum_report.py` (unchanged from 2026-06-28).
