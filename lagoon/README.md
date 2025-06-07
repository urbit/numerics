#   Numerical Libraries for Urbit

**Current Status**

- `/lib/math` is distributed via the [`%yard` desk](https://github.com/urbit/yard).
- We are preparing a release of Lagoon `%real` with SoftBLAS-powered jets for release with Urbit 410K.
- We are preparing a release of Saloon `%real` without jets for release with Urbit 409K.
- We are working on an implementation of [tinygrad](https://tinygrad.org/) for Maroon.

- [ ] some notion of cursor in an array would be really helpful

---

We envision four libraries and associated jet code living in this repository:

- `/lib/math` offers basic special function support for floating-point atoms.
- Lagoon (Linear AlGebra in hOON) offers BLAS-like operations (like NumPy's pure matrix operations).
  - Lagoon `%real`s are slated to ship with [410 K](https://github.com/urbit/UIPs/pull/45).
  - [SoftBLAS](https://github.com/urbit/SoftBLAS) provides a reproducible software-defined floating-point implementation of parts of BLAS and LAPACK suitable for jetting Lagoon.
  - `/lib/fixed` provides operations for fixed-precision operations.
- Saloon (Scientific ALgorithms in hOON) offers transcendental functions (like NumPy's transcendental functions, optimizers, etc.).
- Maroon (MAchine LeaRning in hOON) offers machine learning algorithms, starting with tinygrad.

##  Type System

- `%real` IEEE 754 float (currently supported)
- `%cplx` IEEE 754 float/BLAS packed complex (planned)
- `%uint` unsigned integers (currently supported)
- `%int2` signed twos-complement integers (planned)
- `%sint` signed ZigZag integers (planned)
- `%unum` → subdivide into one of:
  - `%post` posits (32-bit or 16-bit) (planned)
  - `%vald` valids (32-bit or 16-bit) (planned)
- `%fixp` fixed-precision (planned)

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

##  Lagoon 410K `%real` Release

The 410 K release candidate for Lagoon provides `%real`-valued array operations equivalent to `@rh`, `@rs`, `@rd`, and `@rq` operations in vanilla Hoon.  (Other types have been removed for this release but will be re/introduced in a subsequent release.)  The following arms are provided:

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
- `++mmul`
- `++abs`
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
- `++mpow-n`
- `++is-close`
- `++any` (note boolean)
- `++all` (note boolean)
- `++fun-scalar` (helper function)
- `++trans-scalar` (helper function)
- `++el-wise-op` (helper function)
- `++bin-op` (helper function)

The Hoon release is available in `urbit/numerics`, branch `sigilante/reals-only`.  The Hoon release consists of the following files:

- `/sur/lagoon` for data types necessary to use Lagoon.
- `/lib/lagoon` for operations.
- `/tests/lib/lagoon` for various array operation behavior tests.

A PR is at [#6971](https://github.com/urbit/urbit/pull/6971).

The Vere release is available in `urbit/vere`, branch `sigilante/lagoon-jets`.  The Vere release contains jets for 28 arms.  These have been tested for correctness against the reference Hoon results.

A PR is at [#638](https://github.com/urbit/vere/pull/638).

Points for discussion:

1. Currently we ship Lagoon as `/lib/lagoon` and jet into the `f`/`hex` core in `tree.c`.  Since that generally deals with Zuse-level arms, is it advisable to introduce another level `g`/`hep` for `/lib` jets?
2. `/sur/lagoon` still lists other types in `+$kind`, but these are commented out for the time being.

Nonobvious points to note:

1. The comparison gates for Lagoon flip back to boolean rather than loobean results.  Furthermore, they result in numerical ones (e.g. `0x3f80.0000` for `@rs`) rather than simple `0x1`s.  This is because we want sparse matrices to remain sparse when we eventually support them, and because we want multiplication times the result of a logical operation to set or clear fields appropriately without needing to change the `kind`.  (No solution appears to be completely satisfactory.)
2. `++submatrix` and `++stack` are not jetted yet.  These are both dicey jets to get right due to multiple offsets.  Fortunately, once we have them correct they should work for all `kind`s since they only depend on `bloq` size not `kind`.
3. The rounding mode for `%real` may be set for the core using the `++lake` gate.  This returns a copy of the Lagoon `++la` core with rounding mode changed to one of `?(%n %u %d %z)`.
```hoon
> (cumsum:(lake:la %u) (en-ray:(lake:la %u) [~[7 1] 5 %real ~] ~[.1 .5 .-5 .2 .3 .-20 .-1]))
[meta=[shape=~[1 1] bloq=5 kind=%real fxp=~] data=0x1.c170.0000]
```

---

to make:

- [ ] logspace
- [ ] tensordot
- [ ] bitwise ops
- [ ] eq, ne
- [ ] isnan, isinf (±)
- [ ] pad
- [ ] pow, exp, log, whatever not in Saloon

---

THREE POSSIBILITIES:

1. Logical loobean (0 true, terrible for sparse matrices)

2. Logical boolean (1 true, normal but out of step w/ Hoon)

3. Numeric boolean (0x3f80.0000 for 1, etc.)

Here we follow the third option in `u3qf_la_gth_real()` etc.
