#   Numerical Libraries for Urbit

**Current Status**

- We are preparing a release of Lagoon `%real` with SoftBLAS-powered jets for release with Urbit 410K.

---

We envision three libraries and associated jet code living in this repository:

- Lagoon (Linear AlGebra in hOON) will offer BLAS-like operations (like NumPy's pure matrix operations).
  - Lagoon `%real`s are slated to ship with [410 K](https://github.com/urbit/UIPs/pull/45).
  - [SoftBLAS](https://github.com/urbit/SoftBLAS) provides a reproducible software-defined floating-point implementation of parts of BLAS and LAPACK suitable for jetting Lagoon.
  - `/lib/fixed` provides operations for fixed-precision operations.
- Saloon (Scientific ALgorithms in hOON) will offer transcendental functions (like NumPy's transcendental functions, optimizers, etc.).
  - [`/lib/math`](https://github.com/sigilante/libmath) provides a reference interface in Hoon atoms to which Saloon should adhere.
- Maroon (MAchine LeaRning in hOON) will offer machine learning algorithms as a sidecar to Urbit.

---

- `%real` IEEE 754 float
- `%cplx` IEEE 754 float/BLAS packed complex
- `%uint` unsigned integers
- `%int2` signed twos-complement integers
- `%sint` signed ZigZag integers
- `%unum` → subdivide into one of:
  - `%post` posits (32-bit or 16-bit)
  - `%vald` valids (32-bit or 16-bit)
- `%fixp` fixed-precision

---

Fixed-Point Library `/lib/fixed`

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

TODOs:

- %stack
- %cumsum
- %argmin
- %ravel
- %argmax
- %range
- %submatrix

to make:

- [ ] logspace
- [ ] eq, ne

---

THREE POSSIBILITIES:

1. Logical loobean (0 true, terrible for sparse matrices)

2. Logical boolean (1 true, normal but out of step w/ Hoon)

3. Numeric boolean (0x3f80.0000 for 1, etc.)

Here we follow the third option in `u3qf_la_gth_real()` etc.
