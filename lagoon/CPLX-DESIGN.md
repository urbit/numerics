# Design: `%cplx` complex arrays + `/lib/complex`

Status: **design proposal**, 2026-06-01.  Implementation is a follow-up.

This document specifies a complex-number scalar kind (`%cplx`) for Lagoon and
a supporting pure-Hoon `/lib/complex`, with an in-memory atom layout that is
**bit-identical to SoftBLAS / BLAS complex arrays** so the eventual jet is a
zero-conversion cast.  The motivating consumer is eigenvalue / eigenvector
decomposition in **Saloon**: the eigenvalues of a real matrix are in general
complex, so Saloon's `eig` cannot exist until Lagoon can carry complex arrays.

The atom layout and arithmetic in this document were validated on `~hex`
(`(1+2i)·(3+4i) = -5+10i`, `(1+2i)+(3+4i) = 4+6i`, real-low / imag-high `@rs`
components) before writing.

---

## 1. Grounding: what SoftBLAS actually does

From the vendored SoftBLAS (`vere/ext/softblas`) and the Lagoon jet C
(`numerics/lagoon/vere/noun/jets/i/lagoon.c`):

- **Complex scalar types** (`softblas.h`): `complex16_t`, `complex32_t`,
  `complex64_t`, `complex128_t`, each a struct `{ floatN_t real; floatN_t imag; }`
  — **real field first**, then imag, contiguous.
  - `complex16_t` = 2×16b = **32 bits**
  - `complex32_t` = 2×32b = **64 bits**  (single complex)
  - `complex64_t` = 2×64b = **128 bits** (double complex)
  - `complex128_t` = 2×128b = **256 bits**
- **Storage is interleaved**: an array is `[re0, im0, re1, im1, …]`, stride
  measured in whole complex structs (`incX` in units of `complexN_t`).
- **Routines present**: Level-1 complex — `c/i/z/v` prefixes for
  single/half/double/quad — `*axpy`, `*copy`, `*dotc`, `*scal`, `*swap`,
  `*rot`.  Note `*dotc` is the **conjugate** dot (Σ conj(x)·y); there is **no
  unconjugated `*dotu`**.  **There is no complex GEMM** (no `cgemm`/`zgemm`);
  only real `sgemm`/`dgemm`/… exist.
- **Jet marshalling**: the Lagoon jet reads the data atom with
  `u3r_bytes(0, syz, buf, atom)` where `syz = len * 2^(bloq-3)` bytes, appends
  the `0x1` pin for results, and **casts `buf` directly to `(floatN_t*)`**,
  selecting the C type by `bloq` (`4→f16, 5→f32, 6→f64, 7→f128`).  Bytes are
  copied as-is in native (little-endian) order.
- **Saloon** (`numerics/saloon/desk/lib/saloon.hoon`) today is scalar
  elementwise transcendentals/algebra over `%i754` rays.  No eig, no complex.

The consequence that drives the whole design: **if a Lagoon `%cplx` element is
exactly one `complexN_t` (one `2^bloq`-bit slot, real in the low half, imag in
the high half), the jet needs no conversion — it casts the data atom straight
to `complexN_t*` and calls SoftBLAS.**

---

## 2. Atom layout (the core decision)

**`bloq` = log₂ of the *total* complex element width**, i.e. one Lagoon
element = one whole complex number = one SoftBLAS `complexN_t`.  Each complex
element gets its own aura in a new **`@c` (complex) family**, paralleling `@r`
(real): a `@c?` value is a packed pair of the corresponding `@r?` components.

| `%cplx` bloq | element aura | element bits | SoftBLAS type | components | component bloq |
|---|---|---|---|---|---|
| 5 | `@ch` | 32  | `complex16_t`  | 2 × `@rh` (half)   | 4 |
| 6 | `@cs` | 64  | `complex32_t`  | 2 × `@rs` (single) | 5 |
| 7 | `@cd` | 128 | `complex64_t`  | 2 × `@rd` (double) | 6 |
| 8 | `@cq` | 256 | `complex128_t` | 2 × `@rq` (quad)   | 7 |

(Naming follows IEEE component precision: `@cs` = complex-*single* = two 32-bit
floats = a 64-bit element; `@cd` = complex-*double* = two 64-bit floats = a
128-bit element.  These are the BLAS `c` and `z` types respectively.)  Defining
real `@c` auras also lets `get-term` return `%ch/%cs/%cd/%cq` and a future
printer render `re±imi` rather than raw hex.

Within one element slot, **real occupies the low `2^(bloq-1)` bits, imag the
high `2^(bloq-1)` bits**:

```
element  =  (imag << w) | real      where  w = 2^(bloq-1) = (bex (dec bloq))
real     =  (end [0 w] element)
imag     =  (rsh [0 w] element)
```

This matches `complexN_t{real;imag}` over little-endian bytes (real at byte
offset 0 = low bits), so an array of these slots **is** a BLAS interleaved
complex array, and the existing jet `2^(bloq-3)`-byte/element math already
handles the widths (a bloq-6 element = 8 bytes = one `complex32_t`).

Why total-width-bloq rather than component-width-bloq: it keeps Lagoon's "one
element = one `2^bloq` slot" invariant (`cut`/`sew`/`rep`/`get-item` unchanged),
makes the jet a pure cast, and means `dot`/`mmul`/`cumsum` shapes are counted in
complex elements, exactly like every other kind.

**`meta.tail` is unused for `%cplx`** (everything derives from `bloq`): the
component float width is `bloq-1` and the component aura follows from it.  This
is simpler than `%fixp`, which needed `prec` in `tail`.

`zeros` is already correct (all-zero data = complex 0).  `one` (for
`eye`/`ones`) = complex `1.0+0i` = the component float's `1.0` bit pattern in
the low half with zero imag, i.e. just the float-one value (no shift needed).

---

## 3. `/lib/complex` (pure Hoon)

A generic core specialised per width, mirroring `/lib/unum`'s `rpb/rph/rps/rpd`
pattern.  Door names match the aura: `++ch` (`@ch`), `++cs` (`@cs`),
`++cd` (`@cd`), `++cq` (`@cq`); generic `++cx |_ [rnd=rounding-mode]`
parameterised by rounding mode (like the float doors).  Each arm takes/returns
**packed** `@c?` atoms (low=re, high=im) so Lagoon dispatch is a one-liner,
exactly as `add:rpb:unum` is for posits.

Pack/unpack (per the §2 formulas), then component arithmetic via the float door
(`~(add rs rnd)` etc.):

| arm | definition | notes |
|---|---|---|
| `add`  | `[ar+br, ai+bi]` | |
| `sub`  | `[ar-br, ai-bi]` | |
| `mul`  | `[ar·br − ai·bi, ar·bi + ai·br]` | Gauss/Karatsuba 3-mult variant optional |
| `div`  | `[(ar·br+ai·bi)/d, (ai·br−ar·bi)/d]`, `d=br²+bi²` | use **Smith's algorithm** (scale by max(\|br\|,\|bi\|)) to avoid intermediate overflow |
| `neg`  | `[-ar, -ai]` | |
| `conj` | `[ar, -ai]` | |
| `abs`  | `[hypot(ar,ai), 0]` | modulus, returned as a real-valued complex; use `hypot`, not naïve `sqrt(r²+i²)`, to avoid overflow |
| `equ`  | `ar=br ∧ ai=bi` → `?` | |
| `neq`  | `¬equ` | |
| `arg`  | `[atan2(ai,ar), 0]` | optional; needs real `atan2` |

**No total order.**  ℂ is not ordered; `gth`/`gte`/`lth`/`lte` are intentionally
**not defined** and must crash with a clear message
(`~|('complex: no total order; use abs/equ' !!)`).  Only `equ`/`neq` are
total.  (Alternative considered and rejected as default: order-by-modulus —
it silently makes `max`/`min` "work" with a non-standard meaning.  If wanted,
add an explicit `max-abs` later rather than overloading `gth`.)

Decimal `sqt` and `atan2` for `abs`/`arg` come from the same float engine
Saloon's `sqrt` already uses.

---

## 4. Lagoon `%cplx` wiring

Mirror the `%int2`/`%unum`/`%fixp` integration (same six dispatch sites):

- **`sur/lagoon`**: uncomment `%cplx` in the `kind` union.
- **`fun-scalar`** (`=, meta`, dispatch on `bloq` 5/6/7/8 → component door):
  `add`/`sub`/`mul`/`div` via `/lib/complex`; `equ`/`neq` return complex
  `one`/`zero`; **`gth`/`gte`/`lth`/`lte` crash** (no order); no `%mod`/`%pow`.
- **`trans-scalar`**: `%abs` → modulus (as real-valued complex); propose adding
  a new `%conj` op to the `ops` union (other kinds `!!` on `%conj`).
- **`get-term`**: returns the element aura `%ch/%cs/%cd/%cq`.  Until a complex
  printer is registered these render as raw hex, but the aura is now nameable
  (vs `%unum`/`%fixp`, which fall back to `%ux`).
- **`eye`/`ones`**: constant `one` = component-float `1.0` (imag 0).
- **`scale`**: raw `@ux` pack, like `%int2`/`%unum`/`%fixp`.
- **`change`/convert**: guard with `~|` for now.  Eventual map: `%i754`→`%cplx`
  (real r ↦ r+0i), `%cplx`→`%i754` (real part, or `abs`).

### Reductions

Complex float arithmetic has **no exact accumulator** here (unlike the posit
quire or fixed-point integer sums), so `%cplx` reductions round per operation,
like `%i754` — do **not** claim fused exactness.

- **`cumsum`**: independent float sums of the real and imag parts.
- **`dot`**: `Σ aᵢ·bᵢ`, **unconjugated**, to stay consistent with the real
  `dot` across kinds.  Pure-Hoon now.
- **`dotc`** (propose as a new arm): `Σ conj(aᵢ)·bᵢ`, the Hermitian inner
  product Saloon's eig/QR actually wants.  This is the one that maps onto
  SoftBLAS `*dotc` for a future jet (note: there is no `*dotu`, so the
  unconjugated `dot` jet would need pure-Hoon fallback or a `*dotu` added to
  SoftBLAS).
- **`mmul`**: per-cell complex dot, pure-Hoon (like `mmul-unum`).  **No complex
  GEMM exists in SoftBLAS**, so the eventual jet must build it from real GEMMs:
  for `A=Aᵣ+iAᵢ`, `B=Bᵣ+iBᵢ`, either the **4-real-GEMM** form
  (`re = AᵣBᵣ − AᵢBᵢ`, `im = AᵣBᵢ + AᵢBᵣ`) or the **3-mult Karatsuba** form,
  over `sgemm`/`dgemm`, after de-interleaving the complex array into planar
  real/imag buffers and re-interleaving the result.

---

## 5. Jets (later; not this PR)

- **Level-1** (`axpy`/`scal`/`copy`/`swap`/`dotc`): the data atom casts
  directly to `complexN_t*` (§2), so these are thin wrappers over SoftBLAS
  `c/i/z/v` routines — the cheapest jets in the whole numerics tree.
- **`dot` (unconjugated)**: no SoftBLAS primitive; pure-Hoon or add `*dotu`.
- **`mmul`**: 4-real-GEMM or 3M Karatsuba over real `sgemm`/`dgemm` (§4); needs
  de-interleave/re-interleave helpers.
- Consistent with the rest of the tree, `%cplx` (like `%unum`/`%fixp`) runs
  pure-Hoon until a consumer needs the speed.

---

## 6. The payoff: Saloon eigen

`%cplx` is the enabling type for Saloon's spectral layer:

- Real symmetric / Hermitian matrices: real eigenvalues, but eigenvectors and
  intermediate Givens/Householder rotations benefit from complex support.
- General real matrices: complex-conjugate eigenvalue pairs — **require**
  complex arrays for `eigvals`/`eigvecs` at all.

Dependency chain: **`/lib/complex` + Lagoon `%cplx` (this design) → complex
`dotc`/`mmul` → Saloon `eig`** (Hessenberg reduction + shifted QR algorithm).
The QR algorithm is the large follow-on effort; it is unblocked, not done, by
this work.

---

## 7. Decisions

Resolved (2026-06-05):

1. **Auras** — dedicated `@c` family: `@ch`/`@cs`/`@cd`/`@cq` (§2).
2. **Ordering** — `gth`/`gte`/`lth`/`lte` **crash**; only `equ`/`neq` defined.
3. **`%conj` op** — **add** to the shared `ops` union (clean `trans-scalar`
   path); also exposed through `/lib/complex`.
4. **Widths first** — ship `@cs` (bloq 6) and `@cd` (bloq 7), the BLAS `c`/`z`
   types; `@ch`/`@cq` (bloq 5/8) are additive follow-ups.

Still under discussion:

5. **`dot` conjugation** — bilinear `dot` (= numpy `@`/`dot`, consistent with
   `mmul`) plus a separate Hermitian `dotc` (= numpy `vdot`, BLAS `*dotc`), vs
   making `%cplx`'s `dot` itself Hermitian.  See §4.
6. **`change` semantics** — specifically `%cplx`→real: real part (recommended,
   = numpy `.real`) vs modulus (already available as `abs`).  See §4.

---

## 8. Phasing

- **Phase 1** (this design → next PR): `/lib/complex` pure Hoon (add/sub/mul/div/
  neg/conj/abs/equ/neq, oracle-checked) + Lagoon `%cplx` element-wise +
  pure-Hoon `dot`/`dotc`/`mmul`/`cumsum`.  Tests mirror `/tests/lib/lagoon-fixp`.
- **Phase 2**: SoftBLAS jets — Level-1 direct-cast; `mmul` via real-GEMM
  decomposition.
- **Phase 3**: Saloon `eig` (Hessenberg + shifted QR) consuming `%cplx`.

An offline oracle (`tools/complex_check.py`, NumPy `complex64`/`complex128`) is
the ground truth for Phase 1, the same pattern as `posit_check.py` /
`fixed_check.py`.
