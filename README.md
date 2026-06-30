#   Numerical Libraries for Urbit

**Status ~2026.6.27:  Two's complement and unum/posit/quire support added; Lagoon updated for Vere64 by @matthewlevan; Saloon conversion to Chebyshev basis functions completed by @sigilante.  Work proceeds on unum/posit/quire jetting using SoftUnum.**

![An evocative scene of a mysterious futuristic castle in the style of Flash Gordon](./img/hero-scene.jpg)

This repository organizes the core numerical computing apparatus for Urbit:

- `libmath` provides scalar arithmetic libraries; all live in `libmath/desk/lib/`.
  - [`README.md`](./libmath/README.md)
  - `/lib/math` — four-precision transcendentals (`@rs`/`@rd`/`@rh`/`@rq`), jetted via SoftFloat.
  - `/lib/unum` — 2022 Posit Standard (`@rpb`/`@rph`/`@rps`/`@rpd`/`@rpq`) with quire.
  - `/lib/complex` — BLAS-interleaved complex numbers (`@ch`/`@cs`/`@cd`/`@cq`).
  - `/lib/fixed` — fixed-point Q-format arithmetic.
  - `/lib/twoc` — two's-complement signed integers.
- Lagoon (Linear AlGebra in hOON) offers operations in the tradition of BLAS and LAPACK (like NumPy's pure matrix operations).
  - [`README.md`](./lagoon/README.md)
  - `lagoon/desk` contains the Hoon-specific code for Lagoon.
    - `/lib/lagoon` is the main library for Lagoon operations.
    - `/sur/lagoon` supplies type headers for Lagoon.
  - `lagoon/vere` contains the C jets for the Vere runtime.
- Saloon (Scientific ALgorithms in hOON) affords element-wise transcendentals and eigendecomposition over Lagoon rays.
  - [`README.md`](./saloon/README.md)
  - `saloon/desk` contains the Hoon-specific code for Saloon.
    - `/lib/saloon` is the main library for Saloon operations.
- Maroon (MAchine LeaRning in hOON) implements machine learning algorithms as a sidecar to Urbit.
  - [`README.md`](./maroon/README.md)
  - `/desk` contains the Hoon-specific code for Maroon, currently an in-progress tinygrad implementation.

The Urbit Foundation also provides [SoftBLAS](https://github.com/urbit/SoftBLAS) to support software-defined jetting.  It is used in the Lagoon jets.
