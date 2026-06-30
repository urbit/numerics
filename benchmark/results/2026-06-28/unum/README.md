# `/lib/unum` (posit) jet benchmark — 2026-06-28

Three-way per-call comparison of the `/lib/unum` arms: **interpreted** Hoon vs
**Python/SoftUnum** (the C lib via ctypes) vs **jetted** Hoon, at posit8/16/32
(`rpb`/`rph`/`rps`).  Full table in [`table.txt`](table.txt).

## Method (mirrors the `/lib/math` `bench-math` protocol)

- `gen/bench-unum-grid` (→ `lib/unum-cells`) precomputes a posit input list
  *outside* the timed loop, then folds each arm `n` times inside `~>(%bout)`,
  which slogs `took …`.  Per-call = `(arm − base) / n`.  The `%bout` value with
  the dots stripped is microseconds.
- **One jet binary, hint toggle.**  Jetted = `/lib/unum` with its `~%`/`~/`
  hints; interpreted = the same file with the hints commented out (so no jet
  matches and the pure Hoon runs).  Both on a hoon-135 fakezod (the 408k pill).
- Jetted measured at `n=100,000`; interpreted at `n=100` (interpreted posit
  transcendentals are ~20–60 ms/call, so a larger `n` is impractical).
- Python column: `tools/bench_unum_report.py` times SoftUnum via ctypes
  (200k calls/arm) — the raw C speed a Python user sees (call overhead included).

## Headline numbers (jetted vs interpreted)

| class | jetted | interpreted | speedup |
|---|---|---|---|
| arithmetic (`add`/`sub`/`mul`/`div`/`fma`) | ~1.5–2 µs | ~270–530 µs | **~180–295×** |
| `sqrt` posit8 | 1.4 µs | 338 µs | 238× |
| transcendentals (`exp`/`log`/`sin`/`cos`/`pow`) | ~12–60 µs | ~18–57 ms | **~760–2000×** |
| `fdp` (fused dot product, per element) | 0.09–0.27 µs | ~205 µs | **760–2300×** |

The fused dot product shows the largest win: the jet runs the whole quire
accumulation in C (no per-element noun allocation), so posit8 `fdp` is ~0.09
µs/element — ~2300× the interpreted Hoon.

## Caveat: posit16/32 `sqrt`/`atan` are slow even jetted

`sqt`/`atan` widen sharply at posit16/32 because SoftUnum's wide path uses a
512-bit fixed `wide_t` with a bit-by-bit integer sqrt; e.g. jetted `atan:rps` is
~850 µs/call (still 75× the interpreted Hoon, but far off the posit8 45 µs).
The 512-bit `isqt`/AGM is the bottleneck — a candidate for a faster wide sqrt if
posit32 transcendental throughput matters.

## Files

- `table.txt` — the full per-call table.
- `jetted-n100000.txt`, `interp-n100.txt` — raw scraped `%bout` grids.
- harness: `benchmark/desk/lib/unum-cells.hoon`, `benchmark/desk/gen/bench-unum{,-grid}.hoon`,
  `benchmark/tools/bench_unum_report.py`.
