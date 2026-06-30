#   Numerical Libraries for Urbit

**Current Status**

Lagoon ships six element kinds — `%i754` (IEEE 754 floats at @rh/@rs/@rd/@rq),
`%uint` (unsigned integers), `%int2` (two's-complement signed integers via
`/lib/twoc`), `%unum` (2022 Posit Standard via `/lib/unum`), `%cplx`
(BLAS-interleaved complex floats via `/lib/complex`), and `%fixp` (Q-format
fixed-point via `/lib/fixed`).  Array operations are jetted for `%i754` via
SoftBLAS; all other kinds are pure-Hoon.  Saloon provides element-wise
transcendentals and symmetric/Hermitian eigendecomposition (`++eig`) over Lagoon
rays.

---

The numerics repository provides:

- `/lib/math` — scalar transcendentals for `@rs`/`@rd`/`@rh`/`@rq`, jetted via SoftFloat.
- Lagoon — BLAS-like N-D array operations.  [SoftBLAS](https://github.com/urbit/SoftBLAS) provides reproducible software-defined FP for jetting the `%i754` kind.
- Saloon — element-wise transcendentals and eigendecomposition over Lagoon rays.
- Supporting scalar libraries: `/lib/unum`, `/lib/complex`, `/lib/fixed`, `/lib/twoc`.

##  Type System

The element `kind` lives in `/sur/lagoon` (the old `%real` is now `%i754`):

- `%i754` IEEE 754 float — `@rh`/`@rs`/`@rd`/`@rq` (supported)
- `%uint` unsigned integers (supported)
- `%int2` signed two's-complement integers, `/lib/twoc` (supported)
- `%unum` unum/posits — `@rpb`/`@rph`/`@rps`/`@rpd`/`@rpq`, `/lib/unum` (supported); `@rpq` (posit-128, bloq 7) is a non-standard extension beyond the 2022 Posit Standard; `++get-term` does not yet handle bloq 7 for `%unum`
- `%cplx` BLAS-packed complex — `@ch`/`@cs`/`@cd`/`@cq`, `/lib/complex` (supported)
- `%fixp` fixed-point Q a.b, `/lib/fixed`; precision `[a b]` in `meta.tail` (supported)

(`%sint` ZigZag integers and `%vald` valids remain possible future additions.)
Arms added since the `%real`-only release below include `++dotc` (Hermitian dot)
and `++conj` (elementwise conjugate); Saloon adds eigendecomposition (`++eig`).

##  Fixed-Point Library `/lib/fixed`

```
> `@ub`(add:fixed 0b100.0001.0000 [8 8] 0b101.0001.0000 [8 8])
0b1001.0010.0000

> `@ub`(sub:fixed 0b100.0001.0000 [8 8] 0b101.0001.0000 [8 8])
0b1.1111.1111.0000.0000

> `[@ub ^]`(mul:fixed 0b100.0001.0000 [8 8] 0b101.0001.0000 [8 8])
[0b1.0100.1001.0001.0000.0000 17 16]

> `[@ub ^]`(div:fixed 0b100.0001.0000 [8 8] 0b101.0001.0000 [8 8])
[0b0 17 16]

> `[@ub ^]`(div:fixed 0b1.0100.0001.0000 [8 8] 0b101.0001.0000 [8 8])
[0b11 17 16]
```

##  Current Arms (`%i754` and all kinds)

The following arms are provided:

- `++print`
- `++slog`
- `++to-tank`
- `++get-term`
- `++squeeze`
- `++submatrix`
- `++product`
- `++gather`
- `++get-item`
- `++set-item`
- `++get-row`
- `++set-row`
- `++get-col`
- `++set-col`
- `++get-bloq-offset`
- `++get-item-number`
- `++strides`
- `++get-dim`
- `++get-item-index`
- `++ravel`
- `++en-ray`
- `++de-ray`
- `++get-item-baum`
- `++fill`
- `++spac` (helper function)
- `++unspac` (helper function)
- `++scalar-to-ray`
- `++eye`
- `++zeros`
- `++ones`
- `++iota`
- `++magic`
- `++range`
- `++linspace`
- `++urge`
- `++scale`
- `++max`
- `++argmax`
- `++min`
- `++argmin`
- `++cumsum`
- `++prod`
- `++reshape`
- `++stack`
- `++hstack`
- `++vstack`
- `++transpose`
- `++diag`
- `++trace`
- `++dot`
- `++dotc` — Hermitian (conjugate) dot product
- `++mmul`
- `++mmul-unum` — matrix multiply for `%unum` arrays (via quire)
- `++mmul-fixp` — matrix multiply for `%fixp` arrays
- `++abs`
- `++conj` — element-wise conjugate
- `++add-scalar`
- `++sub-scalar`
- `++mul-scalar`
- `++div-scalar`
- `++mod-scalar`
- `++add`
- `++sub`
- `++mul`
- `++div`
- `++mod`
- `++pow-n`
- `++gth` (note boolean)
- `++gte` (note boolean)
- `++lth` (note boolean)
- `++lte` (note boolean)
- `++equ` — element-wise equality (numeric boolean)
- `++neq` — element-wise inequality (numeric boolean)
- `++mpow-n`
- `++is-close`
- `++any` (note boolean)
- `++all` (note boolean)
- `++change` — convert between element kinds
- `++fun-scalar` (helper function)
- `++trans-scalar` (helper function)
- `++el-wise-op` (helper function)
- `++bin-op` (helper function)

Lagoon is shipped in `urbit/urbit` (Hoon: `/lib/lagoon`, `/sur/lagoon`) and `urbit/vere` (C jets via SoftBLAS).  All six element kinds are active in `+$kind`; none are commented out.

Nonobvious points to note:

1. The comparison gates for Lagoon flip back to boolean rather than loobean results.  Furthermore, they result in numerical ones (e.g. `0x3f80.0000` for `@rs`) rather than simple `0x1`s.  This is because we want sparse matrices to remain sparse when we eventually support them, and because we want multiplication times the result of a logical operation to set or clear fields appropriately without needing to change the `kind`.  (No solution appears to be completely satisfactory.)
2. `++submatrix` and `++stack` are not jetted yet.  These are both dicey jets to get right due to multiple offsets.  Fortunately, once we have them correct they should work for all `kind`s since they only depend on `bloq` size not `kind`.
3. The rounding mode for `%i754` may be set for the core using the `++lake` gate.  This returns a copy of the Lagoon `++la` core with rounding mode changed to one of `?(%n %u %d %z)`.
```hoon
> (cumsum:(lake:la %u) (en-ray:(lake:la %u) [~[7 1] 5 %i754 ~] ~[.1 .5 .-5 .2 .3 .-20 .-1]))
[meta=[shape=~[1 1] bloq=5 kind=%i754 fxp=~] data=0x1.c170.0000]
```

---

to make:

- [ ] logspace
- [ ] tensordot
- [ ] bitwise ops
- [ ] isnan, isinf (±)
- [ ] pad
- [ ] pow, exp, log, whatever not in Saloon

