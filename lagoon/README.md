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
- `++sum-dim` — sum along a dimension (dropping it)
- `++prod-dim` — product along a dimension
- `++max-dim` / `++min-dim` — extremum along a dimension
- `++argmax-dim` / `++argmin-dim` — index of the first extremum along a dimension (`%uint` bloq 6)
- `++mean` / `++mean-dim` — arithmetic mean (`%i754`)
- `++var` / `++var-dim` — variance with a `ddof` (`%i754`)
- `++std` / `++std-dim` — standard deviation (`%i754`)
- `++norm` / `++norm-dim` — `%l1`/`%l2`/`%linf`/`%fro` norms (`$norm-ord`)
- `++sort-dim` — sort each slice along a dimension, `%asc` or `%des`
- `++argsort-dim` — the sorting permutation (`%uint` bloq 6, stable)
- `++argtop-dim` — indices of the `k` largest along a dimension
- `++take-dim` — NumPy's `take_along_axis`; pairs with the two above
- `++cdist-sq` / `++pdist-sq` — squared Euclidean distance matrices (`%i754`)
- `++broadcast-to` — expand to a shape by NumPy broadcasting rules
- `++prod-list` / `++drop-dim` / `++dim-parts` / `++dim-slices` / `++slice-flat` /
  `++fold-slice` / `++slice-op` / `++slice-idx` / `++rank-slice` /
  `++norm-slice` / `++i754-sun` / `++i754-sqt` (helper functions)
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
4. The axis-wise arms (`++sum-dim` and friends) fold LEFT TO RIGHT along the axis, seeded with each slice's first element.  The whole-array `++cumsum`, `++max`, and `++min` fold right to left (they use `+reel`), so on an inexact sum `++sum-dim` over a rank-1 ray and `++cumsum` can differ in the last bits.  The left fold is the order a C or Rust kernel walks, which is why it is the one a jet must reproduce; `/tests/lib/lagoon-axis-rounding` pins it with a vector whose two fold orders disagree.
5. `++cdist-sq` sums squared differences directly rather than via the Gram identity `|x|^2 + |y|^2 - 2*A*B^T`.  The identity is much faster (a single `++mmul`) but rounds differently and can produce small negative entries, so it is not a legal jet for this arm.  The same test file pins that too.
6. The index-producing arms (`++argmax-dim`, `++argmin-dim`, `++argsort-dim`, `++argtop-dim`) return `%uint` rays of bloq 6 whatever the input width, and `++take-dim` consumes them.  Reductions drop the reduced dimension (NumPy's `keepdims=False`), and a rank-1 ray reduces to shape `~[1]`; `++argtop-dim` is the exception, replacing the dimension with `k`.

---

to make:

- [ ] logspace
- [ ] tensordot
- [ ] bitwise ops
- [ ] isnan, isinf (±)
- [ ] pad
- [ ] pow, exp, log, whatever not in Saloon

