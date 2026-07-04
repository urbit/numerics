# SoftBLAS vs. native BLAS `gemm` — 2026-07-03

Single-precision (`sgemm`) and double-precision (`dgemm`) matrix-multiply,
SoftBLAS (deterministic, SoftFloat-based) vs. Apple Accelerate (`vecLib`
`cblas_sgemm`/`cblas_dgemm`, the stock hardware BLAS on this machine), both
single-threaded, at matrix sizes 10–1000. Full table in
[`table.txt`](table.txt); raw per-size output in the four `*-raw.txt` files.
This is the data behind the mss.tex "The price of determinism" section.

## Why this run exists, not a citation of the old numbers

`~urbit/SoftBLAS/benchmarking/results.txt` already has an sgemm curve
("1970x slower than Python Numpy; 517x slower than OpenBlas") but it is from
the `ed948c9`/`1d5c733` commits, both dated **2024-08-22** — machine
unknown, almost two years stale relative to this run, and there is no record
of what CPU or thread count produced it. It also compares against **OpenBLAS**,
which is not installed on this machine (`brew list openblas` fails; no
`libopenblas*` found anywhere under `/opt/homebrew` or `/usr/local`). Rather
than install and tune an OpenBLAS build to reproduce a stale comparison,
this run uses **Apple's Accelerate framework** (`cblas_sgemm`/`cblas_dgemm`
via `<Accelerate/Accelerate.h>`), which is the true zero-install "stock BLAS"
on this exact machine (Apple M4 Max, arm64, macOS 26.5.1) — the same host
the paper's other benchmark tables were measured on. No new package was
installed for this run.

Additionally, the existing `benchmarking/benchmark_sgemm.c` harness **no
longer compiles** against current SoftBLAS: `sgemm`/`dgemm` gained a
trailing `rndMode` parameter (`5a033e5`, "Finish rndMode wiring", 2026-05-29)
and the `svec()` helper the old harness called no longer exists in
`softblas.h`. Both had to be fixed/adapted (see "Changes to SoftBLAS" below)
before anything would build.

## Method

- **SoftBLAS**: `benchmark_sgemm.c`/`benchmark_dgemm.c` in
  `~/urbit/SoftBLAS/benchmarking/`, calling `sgemm`/`dgemm` directly against
  `libsoftblas.a` (built from SoftBLAS commit `932d51b`, 2026-06-06) +
  `SoftFloat/build/Darwin-arm64-clang/softfloat.a`. Rounding mode fixed at
  `'n'` (round-nearest-even) for every call.
- **Native**: `benchmark_sgemm_accelerate.c`/`benchmark_dgemm_accelerate.c`
  (new files, same structure/RNG seed as the SoftBLAS harness), calling
  `cblas_sgemm`/`cblas_dgemm` from `-framework Accelerate`.
- Both sides: row-major, no transpose, random `A`/`B`/`C`/`alpha`/`beta` in
  `[-1, 1]`, fixed `srand(0)` seed, `CLOCK_MONOTONIC_RAW` wall time around
  the single `gemm` call only (allocation/fill excluded), mean ± stderr
  reported over repeated calls per size.
- **Single-threaded**: `VECLIB_MAXIMUM_THREADS=1` and `OMP_NUM_THREADS=1` set
  for the Accelerate binaries. This matters — confirmed directly: at n=1000,
  `sgemm` via Accelerate takes 1.35 ms/call with the thread limit set vs.
  0.92 ms/call without it, i.e. Accelerate does auto-multithread at this size
  when left uncapped. SoftBLAS is a plain triple loop with no threading of
  its own, so no cap was needed on that side.
- Matrix sizes 10, 20, 30, 40, 50, 75, 100, 200, 500, 1000 (extended past the
  old harness's 500 cap since Accelerate finishes n≤500 too fast — sub-µs at
  n=30 — to be a meaningful comparison point on its own). Repeat counts are
  scaled down at large n to bound total wall time: 200 reps for n≤50, 100 for
  n=75, 50 for n=100, 20 for n=200, 5 for n=500, 3 for n=1000 — SoftBLAS's
  `sgemm` alone takes ~16 s/call at n=1000, so the whole SoftBLAS sweep (both
  precisions) already runs to a couple of minutes at these rep counts.
- `clang` (Apple clang 21.0.0), `-O2`, arm64, macOS 26.5.1.

## Headline result: the slowdown is NOT a fixed tax — it grows with n

| n | sgemm SoftBLAS | sgemm Accelerate | **sgemm slowdown** | dgemm SoftBLAS | dgemm Accelerate | **dgemm slowdown** |
|---:|---:|---:|---:|---:|---:|---:|
| 10 | 21.9 µs | 7.3 µs | 3.0× | 20.8 µs | 4.7 µs | 4.4× |
| 20 | 160 µs | 9.0 µs | 17.7× | 151 µs | 1.9 µs | 80.0× |
| 30 | 521 µs | 0.79 µs | 661× | 496 µs | 1.4 µs | 354× |
| 40 | 1.20 ms | 1.4 µs | 873× | 1.16 ms | 2.3 µs | 496× |
| 50 | 2.27 ms | 1.9 µs | 1,202× | 2.22 ms | 4.1 µs | 546× |
| 75 | 7.44 ms | 3.2 µs | 2,322× | 7.25 ms | 7.1 µs | 1,020× |
| 100 | 17.3 ms | 5.8 µs | 2,970× | 17.8 ms | 13.9 µs | 1,282× |
| 200 | 133 ms | 22.5 µs | 5,921× | 135 ms | 57.9 µs | 2,342× |
| 500 | 2.02 s | 203 µs | 9,947× | 2.02 s | 675 µs | 3,000× |
| 1000 | 16.05 s | 1.35 ms | **11,853×** | 16.25 s | 4.99 ms | **3,257×** |

At n=10 the gap is only ~3–4×; by n=1000 it is ~11,900× for `sgemm` and
~3,300× for `dgemm`. This is the central finding: **the determinism tax is
not a constant multiplier, it compounds with problem size**, because the two
implementations are on different asymptotic footing in practice even though
both are nominally O(n³) FLOPs:

- **SoftBLAS scales as a clean n³** with no inflection: `sgemm` time
  ratios between consecutive size jumps track the cube of the size ratio
  almost exactly (e.g. 500→1000 is a 2× size jump and a 7.9× time jump,
  cube-of-2 = 8). This is expected — `sgemm.c`/`dgemm.c` are a direct triple
  loop over SoftFloat `f32_add`/`f32_mul` calls with no cache blocking,
  tiling, or vectorization; every scalar multiply-add pays a full software
  floating-point call.
- **Accelerate scales sub-cubically until it doesn't.** From n=100 to n=200
  (2× size), `sgemm` time only grows 3.9× (vs. 8× for a clean cube) — cache
  blocking and vectorization are absorbing the extra work. From n=500 to
  n=1000, growth is back up to 6.7×, closer to cubic, as the working set
  outgrows the blocking the implementation is tuned for. SoftBLAS has no
  such regime change because it has no blocking to run out of.
- The implied throughput makes the mechanism concrete: at n=1000, SoftBLAS
  sustains **~0.12 GFLOP/s** for both `sgemm` and `dgemm` (precision-
  independent, because SoftFloat call overhead — not memory bandwidth —
  dominates), while Accelerate hits **~1,477 GFLOP/s** for `sgemm` and
  **~401 GFLOP/s** for `dgemm`, single-threaded. The three-orders-of-
  magnitude-plus Accelerate numbers are consistent with Apple Silicon's AMX
  matrix coprocessor, which Accelerate/vecLib dispatches to for `gemm`-shaped
  workloads even from a single CPU thread (AMX is a per-core coprocessor,
  not a multithreading effect) — that's the specific hardware feature the
  "price of determinism" is being paid against.

## Caveat: the smallest sizes are noise-dominated, not signal

At n=10/20, Accelerate's own reported stderr is comparable to its mean
(e.g. n=20 `sgemm`: 9.0 µs ± 8.3 µs) — these calls are fast enough that
timer resolution, first-call/cache warm-up, and scheduling jitter dominate
the measurement, not the actual FLOP cost. Notice `sgemm` accelerate
actually reports a *smaller* mean at n=30 (0.79 µs) than at n=20 (9.0 µs) —
that inversion is measurement noise, not a real effect. Treat the slowdown
factors at n≤20 as order-of-magnitude only; the clean, monotonically
growing trend from n=30 upward (661× → 11,853× for `sgemm`) is the reliable
part of the story and is what the paper should lead with.

## Changes made to SoftBLAS to make this run possible

Not committed (per instructions) — left as a working-tree diff in
`~urbit/SoftBLAS` for review:

- `benchmarking/benchmark_sgemm.c`: fixed to compile against the current
  `sgemm` signature (added the trailing `'n'` rndMode argument) and replaced
  the no-longer-existing `svec()` helper with plain `malloc`+`memcpy` (valid
  because `float32_t` is a struct wrapping a single `uint32_t` with the same
  bit layout as IEEE binary32, so a flat byte copy from a native `float`
  array reinterprets in place). Extended the size list from `{10..500}` to
  `{10..1000}` and reduced the large-n repeat counts so the run stays
  reasonably bounded in time.
- `benchmarking/benchmark_dgemm.c` (new): the same fix/structure, for
  `dgemm`/`float64_t`/`double` — no dgemm benchmark previously existed.
- `benchmarking/benchmark_sgemm_accelerate.c`,
  `benchmark_dgemm_accelerate.c` (new): Accelerate-framework analogs of the
  existing (but unbuildable-here, since OpenBLAS isn't installed)
  `benchmark_sgemm_openblas.c`, same RNG/size/timing structure.

No SoftBLAS library or test code was touched — only the benchmarking
harness. Nothing was committed in either repo.

## Files

- `table.txt` — full sgemm+dgemm table (raw numbers, this README's table
  rounds for readability).
- `sgemm-softblas-raw.txt`, `sgemm-accelerate-raw.txt`,
  `dgemm-softblas-raw.txt`, `dgemm-accelerate-raw.txt` — raw stdout from each
  of the four benchmark binaries (mean ± stderr per size, full precision).
- Harness (in `~/urbit/SoftBLAS/benchmarking/`, uncommitted):
  `benchmark_sgemm.c`, `benchmark_dgemm.c`, `benchmark_sgemm_accelerate.c`,
  `benchmark_dgemm_accelerate.c`.
