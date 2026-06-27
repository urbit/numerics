# Math transcendentals benchmark — 2026-06-26 (32-bit)

Per-call wall-clock (ns), baseline-subtracted, mean over sets 1–4 (set 0 dropped).
Binary: `urbit-jet-32` / `urbit-nojet-32`, both vere 4.4 (`b591168a6f`), on a fresh
hoon-136 fakezod booted from the downloaded `urbit-v4.4.pill`.

## Methodology
- **jetted**: `urbit-jet-32` (math C jet attaches), n=100k.
- **cheb-interp**: `urbit-nojet-32` (same math.hoon, NO math jet registered → clean
  interpreted Chebyshev; jetted SoftFloat base ops), n=10k.  NB: do NOT interpret
  by commenting `~%` on the jet binary — that makes the runtime re-attempt the
  unparentable jet every call (`fund: parent not found` flood + overhead).
- **taylor**: legacy iterative `Σxⁱ/i!` lib (bare door, never jetted), n=10 (it's
  ~12 min/cell at 10k for slow-converging asin/acos).
- Inputs precomputed OUTSIDE `~>(%bout)` so the slow interpreted `sun:rd` isn't timed.

## Headline: @rd (f64) jet speedup
| arm | jetted | cheb-interp | taylor | cheb× | taylor× |
|-----|-------:|------------:|-------:|------:|--------:|
| exp | 1.67µs | 296µs | 11.0ms | 177× | 6,565× |
| log | 1.76µs | 164µs | 13.3ms | 94× | 7,600× |
| sin | 1.95µs | 278µs | 0.94ms | 143× | 480× |
| cos | 1.97µs | 275µs | 0.94ms | 140× | 476× |
| tan | 4.53µs | 241µs | 1.30ms | 53× | 288× |
| atan | 2.23µs | 47.6µs | 109ms | 21× | 48,784× |
| asin | 2.07µs | 54.8µs | 205ms | 27× | 99,130× |
| acos | 2.10µs | 51.2µs | 5.17ms | 24× | 2,459× |
| sqt | 1.92µs | 10.3µs | 73.5µs | 5× | 38× |
| cbt | 2.44µs | 467µs | 111ms | 192× | 45,723× |
| pow | 2.87µs | 656µs | 1.68ms | 229× | 585× |
| pow-n | 2.62µs | 217µs | 219µs | 83× | 83× |
| log-2 | 2.34µs | 164µs | 13.3ms | 70× | 5,668× |
| log-10 | 2.38µs | 165µs | 13.4ms | 69× | 5,638× |
| atan2 | 2.72µs | 56.2µs | 11.5ms | 21× | 4,237× |

@rs (f32) is consistent (e.g. exp 185×/4,574×, asin 17×/27,936×, sqt ~1×/35×).

## Exotic precisions (@rh/@rq): interpreted is infeasible
`@rh`/`@rq` base kernel ops are NOT jetted, so interpreted transcendentals are
~10⁶× slower — a full grid is impossible.  Single salvaged point:
- **taylor `@rh exp` = 70.8 s/call** vs jetted ~1.9µs → **~37,000,000×**.
Jetted covers all 4 doors at ~2µs/call regardless; interpreted @rh/@rq are
captured as single points only (see `bench-interp-exotic.csv`).

## Files
- `bench-timing.csv` — 640 rows: jetted (4 doors) + cheb/taylor (@rd/@rs).
- `bench-summary.tsv` — mean±stddev per cell.
- `bench-interp-exotic.csv` — salvaged interpreted @rh points.
- `raw-32-*.txt` — raw dojo scrollback per mode (re-parse with `tools/parse_raw.py`).

## TODO (not yet run)
- 64-bit word-size axis (`urbit-jet-64`/`urbit-nojet-64`).
- More @rh/@rq interpreted single-points (cheb @rh exp; @rq exp).
- Accuracy pass (`gen/bench-acc.hoon`, ULP cheb-vs-taylor; ground truth via rq_check.c).
