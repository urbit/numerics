# numerics benchmark suite

Reusable harness for measuring **speed** (per-call wall-clock) and **accuracy**
(ULP vs reference) of the numerics primitives across implementations, precisions,
and runtime word sizes. Built originally to compare the math transcendentals
(Taylor vs Chebyshev vs jet); structured so new versions and libraries plug in.

## Layout
```
benchmark/
  MANIFEST.md            provenance of every lib version + how to add one
  README.md              this file
  desk/                  mountable Urbit content (mount into a bench ship's %base)
    lib/                 the libraries under test
      math.hoon          current Chebyshev (jetted + cheb-interp)
      math-taylor.hoon   legacy iterative Taylor (ca42387^)
      unum.hoon fixed.hoon complex.hoon twoc.hoon   pure-Hoon (no jets yet)
      bench-core.hoon    shared %bout timing loop + baseline + input helpers
      bench-domains.hoon per-(lib,arm) input-domain table
    gen/
      bench-math.hoon    one timing cell: (impl,door,arm,n,sets) -> result row
      bench-acc.hoon     one accuracy cell: (impl,door,arm) -> max/mean ULP
      bench-<lib>.hoon    (unum/fixed/complex — same pattern, per-lib input gen)
  tools/
    bench_run.sh         host driver: boots ships in tmux, sends cells, scrapes CSV
    bench_summarize.py   aggregates 5 sets -> mean +/- stddev
  results/               dated run records (CSV + summaries + notes); see results/README.md
```

## Axes (math)
- **impl (3):** taylor · cheb-interp · jetted
- **width (4):** @rh f16 · @rs f32 · @rd f64 · @rq f128
- **word size (2):** 32-bit (`~/urbit/vere`) · 64-bit (`~/urbit/vere-ml64`, `-Dvere64=true`)
- **arm (15):** exp log sin cos tan atan atan2 asin acos sqt cbt pow pow-n log-2 log-10
- **reps:** jetted 100k×5; interpreted 10k×5 (drop set 0 as warm-up)

## How to run
1. **Build four binaries** — `{vere, vere-ml64} × {jet, nojet}`. The jet builds are
   the normal runtimes (full math block in 135/136/137 `tree.c`). For each **nojet**
   build, comment out the single `{ "math", 7, 0, _13X_non__math_d, no_hashes }` line
   in the **136** tree (ships boot hoon-136) and rebuild (per `libmath/vere64/README.md`
   / the build memory). Base `@rX` ops stay jetted; only the transcendental falls back
   to Hoon. Verify each: `~>(%bout (exp:rd:math .~1))` ≈ 1 µs (jet) vs ≈ 250 µs (nojet).
2. **Boot fresh fakezods** (one per binary) in tmux from a kelvin-136 pill, mount
   `%base`, and sync `desk/lib` + `desk/gen` in. Never reuse old piers; never `C-c`
   the serf (drive with `C-u`).
3. **Smoke test** (`tools/bench_run.sh --smoke`): n=1000, sets=2, arms {exp,sin,sqt,
   atan2}, doors {rd,rh}, impls {cheb,taylor}. Confirms the jet toggle (~250× split),
   no memoization (baseline-subtracted time positive/stable), no @rh stall, Taylor
   guard (`fail` not crash on out-of-domain).
4. **Full run** (`tools/bench_run.sh`): writes `results/<date>/bench-timing.csv` and
   `bench-accuracy.csv`; `bench_summarize.py` produces `bench-summary.tsv`.

## Adding a version or library (for future trials)
- **New math poly/jet version:** drop the lib as `desk/lib/<name>.hoon` with a DISTINCT
  top-level shape (so the jet matches only the intended one), add a MANIFEST row, and add
  it to the `impl` switch in `gen/bench-math.hoon`.
- **New library (e.g. a future jetted posit):** copy its lib into `desk/lib/`, write
  `gen/bench-<lib>.hoon` reusing `lib/bench-core.hoon` (the timing loop) with a per-lib
  input generator + arm dispatch, and add an entry to MANIFEST + this README.

## Methodology guards (baked into `bench-core.hoon`)
- **Defeat memoization:** inputs vary every iteration (counter-driven), and each call's
  result is folded into a returned accumulator — neither input nor whole loop is a
  cacheable constant. A memoized cell shows ~0 net time (visible failure).
- **Baseline subtraction:** an identical loop minus the primitive isolates per-call cost.
- **@rh stall:** never accumulate a float (sub-epsilon steps freeze); map a `@ud` counter
  into the domain so f16 samples stay distinct.
- **Non-robust legacy (Taylor):** safe-domain inputs + `mule` guard + iteration cap;
  record `fail`, never abort the run.
