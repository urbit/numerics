# SoftBLAS agent session — complex BLAS for Lagoon `%cplx` (Phase 2)

This is a self-contained brief for a session working on **SoftBLAS** (C, built
with Zig).  Paste it as the opening prompt.  Goal: add the complex primitives
Lagoon's `%cplx` array kind needs, so the Lagoon jets can call them instead of
running pure Hoon.

---

## Mission

Add to our SoftBLAS fork:

1. **`cdotu` / `zdotu`** — *unconjugated* complex dot product `Σ xᵢ·yᵢ`
   (single / double). Only the conjugating `cdotc`/`zdotc` exist today.
2. **`cgemm` / `zgemm`** — complex matrix multiply `C ← α·op(A)·op(B) + β·C`
   (single / double). No complex GEMM exists today; only real
   `sgemm`/`dgemm`/`hgemm`/`qgemm`.

Optional / stretch (do only if cheap and time allows): the half/quad analogues
(`idotu`/`vdotu`, `igemm`/`vgemm`) and complex GEMV (`cgemv`/`zgemv`).

These unblock the Lagoon jet for `%cplx` `dot`/`dotc`/`mmul` (a later,
separate task — see "Downstream", you do NOT need to touch Lagoon here).

---

## Where everything is

- **SoftBLAS source (your worktree):** the `urbit/SoftBLAS` repo. Vere only
  *vendors* it — `vere/ext/softblas/` is a thin Zig build wrapper that fetches a
  tarball pinned in `build.zig.zon`:
  `https://github.com/urbit/SoftBLAS/archive/4e702725445ca24d1fda5785b89c607aeb23e3eb.tar.gz`.
  Clone/checkout `urbit/SoftBLAS` at (or after) commit `4e70272` and work there.
- **The public header** (already in the vendored build output, read-only ref):
  `vere/ext/softblas/zig-out/include/softblas.h`. It is the source of truth for
  existing signatures, the `complexN_t` structs, the `cNN_add/sub/mul/div/conj`
  macros, and the `SB_*_ONE/_ZERO/_I` constants. Mirror its conventions exactly.
- **Behavioral oracle (authoritative for results):** the numerics repo's pure-Hoon
  complex lib, `numerics/libmath/desk/lib/complex.hoon`, and its Lagoon wiring
  `numerics/lagoon/desk/lib/lagoon.hoon` (`%cplx`). Plus the design doc
  `numerics/lagoon/CPLX-DESIGN.md` and the NumPy oracle
  `numerics/libmath/tools/complex_check.py`. Your C results MUST match these
  bit-for-bit (same rounding).
- **Stateless rounding:** every SoftBLAS op already takes `const uint_fast8_t
  rndMode` (a SoftFloat rounding mode) per call. Keep it that way — no globals.
  (This matches the stateless-SoftBLAS migration in vere PR #1021,
  branch `neal/lagoon-jet-fixes`; coordinate, don't regress it.)

---

## Layout contract (must match Lagoon exactly)

A complex value is one struct, real field first:

```c
typedef struct { float32_t real; float32_t imag; } complex32_t;  // 8 bytes
typedef struct { float64_t real; float64_t imag; } complex64_t;  // 16 bytes
```

Arrays are **interleaved**: `[re0, im0, re1, im1, …]`, stride `incX` in whole
`complexNN_t` units. This is identical to Lagoon's packed `@cs`/`@cd` element
(real in the low half, imag in the high half, little-endian), so the Lagoon jet
will `memcpy`/cast the noun's data straight to `complexNN_t*` with no
conversion. Do not change the field order or introduce planar storage.

---

## Deliverable 1 — `cdotu` / `zdotu`

Mirror the existing `cdotc`/`zdotc` exactly, but **without conjugating** the
first vector:

```c
complex32_t cdotu(uint64_t N, const complex32_t *CX, int64_t incX,
                  const complex32_t *CY, int64_t incY, const uint_fast8_t rndMode);
complex64_t zdotu(uint64_t N, const complex64_t *CX, int64_t incX,
                  const complex64_t *CY, int64_t incY, const uint_fast8_t rndMode);
```

`cdotc` does `Σ conj(xᵢ)·yᵢ`; `cdotu` does `Σ xᵢ·yᵢ`. Find `cdotc`'s
implementation in the SoftBLAS source and copy it, dropping the `c32_conj` on
`x`. Accumulate term-by-term with `c32_mul` then `c32_add` (per-op rounding via
`rndMode`) — the Hoon oracle has no wide accumulator, so straight per-op
accumulation is what matches. Negative/zero strides: match `cdotc`'s handling.

---

## Deliverable 2 — `cgemm` / `zgemm`

Full BLAS-3 signature, mirroring `sgemm` but with complex types:

```c
void cgemm(const char transA, const char transB,
           const uint64_t M, const uint64_t N, const uint64_t P,
           const complex32_t alpha, const complex32_t *A, const uint64_t lda,
           const complex32_t *B, const uint64_t ldb,
           const complex32_t beta, complex32_t *C, const uint64_t ldc,
           const uint_fast8_t rndMode);
void zgemm(/* same shape with complex64_t */);
```

`C ← α·op(A)·op(B) + β·C`, where `op(X)` is `X`, `Xᵀ`, or `Xᴴ` per `transA`/
`transB` (`'N'`/`'T'`/`'C'` — note complex adds the conjugate-transpose `'C'`,
which `sgemm` lacks; handle it with `c32_conj`). Match `sgemm`'s layout
convention (row/col-major, `lda`/`ldb`/`ldc` leading dimensions) — read
`sgemm` and follow it precisely.

Two valid implementation strategies:

- **(a) Direct triple loop** (recommended first): accumulate each `C[i,j]` with
  `c32_mul`/`c32_add` and the per-call `rndMode`. Simple, obviously correct,
  and it matches the Hoon oracle's per-op rounding. Ship this first.
- **(b) Real-GEMM decomposition** (optional optimization): split into real/imag
  planes and use `sgemm`/`dgemm` — `Cᵣ = AᵣBᵣ − AᵢBᵢ`, `Cᵢ = AᵣBᵢ + AᵢBᵣ`
  (4 real GEMMs) or 3-mult Karatsuba. Faster but changes rounding/accumulation
  order, so it will NOT be bit-identical to the Hoon oracle — only pursue if you
  also relax the oracle expectation and document it. Default to (a).

The Lagoon consumer always calls with `α=1+0i`, `β=0+0i`, `transA=transB='N'`,
but implement the general signature.

---

## Rounding modes

`rndMode` is a SoftFloat rounding mode (the same `uint_fast8_t` every existing
op takes). Lagoon uses four, mapped by the jet's `_la_rnd`:

| Lagoon | IEEE | SoftFloat |
|--------|------|-----------|
| `%n` | nearest, ties-even | `softfloat_round_near_even` |
| `%u` | toward +∞ | `softfloat_round_max` |
| `%d` | toward −∞ | `softfloat_round_min` |
| `%z` | toward zero | `softfloat_round_minMag` |

Use `rndMode` for every component float op; never read a global.

---

## Verification

1. **SoftBLAS's own test harness:** find it in the SoftBLAS repo (there are
   existing tests for `sgemm`, `cdotc`, etc.). Add cases for the new functions
   in the same style and make them pass.
2. **Cross-check against the numerics oracle.** Canonical vectors (float32 /
   `complex32_t`); the Hoon packed-`@cs` hex (real-low/imag-high) is in
   parentheses so you can compare byte patterns:
   - `c32_mul(1+2i, 3+4i) = -5+10i`  (`0x4120.0000.c0a0.0000`)
   - `cdotu([1+2i, 3+4i], [5+6i, 7+8i]) = -18+68i`  (`0x4288.0000.c190.0000`)
   - `cdotc([1+2i, 3+4i], [5+6i, 7+8i]) =  70-8i`   (`0xc100.0000.428c.0000`)  ← sanity vs existing
   - `cgemm([[1+1i, 2+0i],[0, 1+0i]], [[1+0i, 0],[0, 1+1i]]) = [[1+1i, 2+2i],[0, 1+1i]]`
   Run `python3 numerics/libmath/tools/complex_check.py` to regenerate / extend
   these, and (if you have a ship) check against `/tests/lib/lagoon-cplx` and
   `/tests/lib/complex` in numerics, which already pass with the pure-Hoon path.
3. **float32 single-element components:** `1.0=0x3f800000`, `2.0=0x40000000`,
   `3.0=0x40400000`, `4.0=0x40800000`, `5.0=0x40a00000`, `-5.0=0xc0a00000`,
   `10.0=0x41200000`. A `complex32_t` `re+imi` is `{re_bits, im_bits}`.

---

## Downstream (context only — NOT your task here)

The Lagoon jet (`numerics/lagoon/vere/noun/jets/i/lagoon.c`) will later call
these for `%cplx`: element-wise add/sub via `caxpy`+`cscal`, `mul`/`div` per
element, `abs` via `scnrm2`/`c*_div`, `dot`→`cdotu`, `dotc`→`cdotc`,
`mmul`→`cgemm`. The jet marshals the ray's data atom (`u3r_bytes`,
`syz = len·2^(bloq-3)` bytes, strip the `0x1` pin, cast to `complexNN_t*`),
selecting by `bloq` (`%cplx` bloq 6 → `complex32_t`, bloq 7 → `complex64_t`).
You don't wire the jet; you just make the C functions exist and be correct.

---

## Constraints & done

- Branch in `urbit/SoftBLAS`; keep the **stateless** (rnd-per-call) API; reuse
  the `cNN_*` macros and SoftFloat ops; don't break real BLAS or the existing
  complex Level-1.
- Export the new functions in `softblas.h` following the existing prose/order.
- **Done when:** `cdotu`/`zdotu` and `cgemm`/`zgemm` build under the Zig build,
  pass new SoftBLAS tests, and reproduce the canonical vectors above
  bit-for-bit. Then (handoff) the SoftBLAS hash gets bumped in
  `vere/ext/softblas/build.zig.zon` and the Lagoon `%cplx` jet wired — coordinate
  with vere PR #1021 rather than opening a competing vere change.
