# Benchmark version manifest

Provenance of each math-library implementation under `desk/lib/`. To add a new
version for a future trial: drop the lib in `desk/lib/<name>.hoon`, give it a
distinct top-level structure (so the jet only matches the intended one), and add
a row here.

| version | file | source | algorithm | jetted? |
|---------|------|--------|-----------|---------|
| **Taylor** (legacy) | `desk/lib/math-taylor.hoon` | numerics `ca42387^` (`f764eb3`), pre-Chebyshev | iterative `Σ xⁱ/i!` to `rtol`, no range reduction; `atan` AGM-style; composites (pow/cbt/log-2/log-10/atan2) over the iterative kernels | no — bare `\|%` doors, no `~%`/`~/`, never matches the jet |
| **Chebyshev** (current) | `desk/lib/math.hoon` | numerics `7dad431:libmath/desk/lib/math.hoon` | Cody-Waite reduction + minimax Horner (fdlibm sin/cos/asin kernels); ≤1 ULP faithful | yes — `~% %non` / `~/ %math` hints; the C jet (`libmath/.../jets/i/math.c`) is bit-exact to it |

## Pure-Hoon libraries (no jets yet — perf baselines / jetting candidates)
These are included so the suite can measure their interpreted cost and help
prioritize which to jet next. All are pure Hoon at HEAD (`7dad431`).

| library | file | depends on | call shape (for the harness) |
|---------|------|-----------|------------------------------|
| **unum** (posits) | `desk/lib/unum.hoon` | `/+ twoc` | posit door `pp`; ops over encoded posit atoms `p=@`; quire / `fdp` |
| **fixed** | `desk/lib/fixed.hoon` | `/+ twoc` | `prec`-keyed `++ add/sub/mul`; operands `a/b=@` |
| **complex** | `desk/lib/complex.hoon` | (standalone) | door `cs` (rounding-mode keyed); ops over interleaved `[re im]` atoms `[p=@ q=@]` |
| **twoc** | `desk/lib/twoc.hoon` | (standalone) | two's-complement helper used by unum/fixed |

Each gets its own `gen/bench-<lib>.hoon` reusing the shared timing core
(`lib/bench-core.hoon`): same `~>(%bout ...)` loop + baseline subtraction + varying
counter-driven input, with a per-lib input generator and arm dispatch. Math is the
fully-built reference; the pure-Hoon libs follow the same pattern.

## Notes
- The **jetted** measurement and the **Chebyshev-interpreted** measurement use the
  SAME `math.hoon`; the only difference is whether the runtime carries the math jet
  (jet binary vs no-jet binary — see `README.md`).
- `desk/lib/math.hoon` must stay the exact copy the C jet was built against (same
  battery → jet attaches). When the canonical `libmath/desk/lib/math.hoon` changes,
  re-copy it here and bump the source hash above, OR the jet will fall back to Hoon.
- The Taylor version is **non-robust** (crashes/hangs outside safe domains — see the
  domain table in `desk/lib/bench-domains.hoon`); the harness guards each set with
  `mule` and records `fail` rather than aborting.
