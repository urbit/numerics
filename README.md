#   Numerical Libraries for Urbit

**Status ~2024.5.14:  Lagoon releases with IEEE 754 reals in 410K; Saloon releases with corresponding functionality in 409K; Maroon is in progress for tinygrad.**

![An evocative scene of a mysterious futuristic castle in the style of Flash Gordon](./img/hero-scene.jpg)

This repository organizes the core numerical computing apparatus for Urbit:

- `/lib/math` provides basic single-atom transcendental functions; it supersedes [`sigilante/libmath`](https://github.com/sigilante/libmath).
  - [`README.md`](./libmath/README.md)
- Lagoon (Linear AlGebra in hOON) offers operations in the tradition of BLAS and LAPACK (like NumPy's pure matrix operations).
  - [`README.md`](./lagoon/README.md)
  - `/desk` contains the Hoon-specific code for Lagoon.
    - `/lib/lagoon` is the main library for Lagoon operations.
    - `/lib/twoc` supports two's-complement signed integers.
    - `/lib/fixed` supports fixed-precision operations.
    - `/sur/lagoon` supplies type headers for Lagoon.
  - `/vere` contains the C jets for the Vere runtime.
- Saloon (Scientific ALgorithms in hOON) affords transcendental functions (like NumPy's transcendental functions, optimizers, etc.).
  - [`README.md`](./saloon/README.md)
  - `/desk` contains the Hoon-specific code for Saloon.
- Maroon (MAchine LeaRning in hOON) implements machine learning algorithms as a sidecar to Urbit.
  - [`README.md`](./maroon/README.md)
  - `/desk` contains the Hoon-specific code for Maroon, currently an in-progress tinygrad implementation.

The Urbit Foundation also provides [SoftBLAS](https://github.com/urbit/SoftBLAS) to support software-defined jetting.  It is used in the Lagoon jets.
