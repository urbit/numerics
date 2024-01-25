#   Numerical Libraries for Urbit

We envision three libraries and associated jet code living in this repository:

- Lagoon (Linear AlGebra in hOON) will offer BLAS-like operations (like NumPy's pure matrix operations).
  - [SoftBLAS](https://github.com/urbit/SoftBLAS) provides a reproducible software-defined floating-point implementation of parts of BLAS and LAPACK suitable for jetting Lagoon.
  - `/lib/fixed` provides operations for fixed-precision operations.
- Saloon (Scientific ALgorithms in hOON) will offer transcendental functions (like NumPy's transcendental functions, optimizers, etc.).
  - [`/lib/math`](https://github.com/sigilante/libmath) provides a reference interface in Hoon atoms to which Saloon should adhere.
- Maroon (MAchine LeaRning in hOON) will offer machine learning algorithms as a sidecar to Urbit.

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
