#   Numerical Libraries for Urbit

**Status ~2024.5.1:  Development work is taking place preparatory for Lagoon release with 410K and Saloon release with 409K.**

![An evocative scene of a mysterious futuristic castle in the style of Flash Gordon](./img/hero-scene.jpg)

This repository organizes the core numerical computing apparatus for Urbit:

- `/lib/math` provides basic single-atom transcendental functions; it mirrors [`sigilante/libmath`](https://github.com/sigilante/libmath) which is the canonical version.
- Lagoon (Linear AlGebra in hOON) offers operations in the tradition of BLAS and LAPACK (like NumPy's pure matrix operations).
  - `/desk` contains the Hoon-specific code for Lagoon.
  - `/vere` contains the C jets for the Vere runtime.
- Saloon (Scientific ALgorithms in hOON) affords transcendental functions (like NumPy's transcendental functions, optimizers, etc.).
  - `/desk` contains the Hoon-specific code for Saloon.
- Maroon (MAchine LeaRning in hOON) implements machine learning algorithms as a sidecar to Urbit.
  - `/desk` contains the Hoon-specific code for Maroon, currently an in-progress tinygrad implementation.

The Urbit Foundation also provides [SoftBLAS](https://github.com/urbit/SoftBLAS) to support software-defined jetting.  It is used in the Lagoon jets.
