#   Scientific ALgorithms in hOON

Transcendental and algebraic functions for use with Lagoon `+$ray`s.

We support the following functions and special functions:

- `++add`, $+$ addition (pass-through from Lagoon)
- `++sub`, $-$ subtraction (pass-through from Lagoon)
- `++mul`, $\times$ multiplication (pass-through from Lagoon)
- `++div`, $/$ division (pass-through from Lagoon)
- `++fma`, $\text{fma}$ fused multiply-add
- `++neg`, $-$ unary negation
- `++factorial`, $!$ factorial
- `++abs`, $\text{abs}$ (pass-through from Lagoon)
- `++exp`, $\exp$
- `++sin`, $\sin$
- `++cos`, $\cos$
- `++tan`, $\tan$
- `++pow-n`, $\text{pow}$ to integer power
- `++log`, $\log$ (natural logarithm)
- `++log-10`, $\log_{10}$ (log base-10)
- `++log-2`, $\log_{2}$ (log base-2)
- `++pow`, $\text{pow}$
- `++sqrt`, $\sqrt$ (also `++sqt`)
- `++cbrt`, $\cbrt$ (also `++cbt`)

Logical functions:

- `++lth`, $<$ (pass-through from Lagoon)
- `++lte` $\leq$ (also `++leq`) (pass-through from Lagoon)
- `++gth`, $>$ (pass-through from Lagoon)
- `++gte`, $\geq$ (also `++geq`) (pass-through from Lagoon)
- `++equ`, $=$ (pass-through from Lagoon)
- `++neq`, $\neq$ (pass-through from Lagoon)
- `is-close` (pass-through from Lagoon)
- `all-close` (pass-through from Lagoon `++all`)
- `any-close` (pass-through from Lagoon `++any`)

Linear algebra (over Lagoon arrays):

- `++eig`, eigendecomposition of a **symmetric** (`%i754`) or **Hermitian**
  (`%cplx`) matrix via cyclic Jacobi → `[vals=ray vecs=ray]` (real eigenvalues,
  orthonormal/unitary eigenvectors as columns).  Dispatches on `kind`.
- `++eigvals`, eigenvalues only (1-D ray).
- `++eigvecs`, eigenvectors only.

Solvers and factorizations (`%i754` only, bloq 4/5/6/7):

- `++chol`, the Cholesky factor `L` of a symmetric positive-definite matrix
  (`A = L*L^T`, lower triangular); `++chol-unit` is the same thing returning
  `(unit ray)`, with `~` instead of a crash when the matrix is not positive
  definite.
- `++trsv-lo` / `++trsv-up`, forward and back substitution against `L` (the
  latter reads `L` transposed rather than materializing `L^T`).
- `++chol-solve`, `A*x = v` by one factorization plus the two substitutions.
- `++qr`, the thin Householder QR `A = Q*R` (rows >= cols): `Q` has orthonormal
  columns, `R` is upper triangular with an exactly-zero subdiagonal. Signs
  follow LAPACK (`alpha = -sign(x0)*|x|`), so they match `numpy.linalg.qr`,
  negative diagonal entries included. `++trsv-r` is its back substitution.
- `++lstsq`, least squares `min |A*x - v|` as `R^-1 * Q^T*v`, avoiding the
  squared condition number of the normal equations. Needs full column rank.
- `++gram` (`x^T*x`) and `++matvec-t` (`m^T*x`), both without materializing a
  transpose — and so without the lagoon transpose jet, which crashes on
  runtimes older than urbit/vere#1057.
- `++cg`, conjugate gradient from `x0 = 0` → `[x iter rnorm]`, so the caller can
  tell convergence from exhaustion; `++pcg` adds the Jacobi (diagonal)
  preconditioner.  Each iteration is one `++matvec` and two `++dotv`, touching
  the matrix only through products.
- `++svd`, the thin singular value decomposition by **one-sided Jacobi** →
  `[u s v]` with `A = U*diag(S)*V^T` and `s` sorted DESCENDING (unlike `++eig`,
  whose order is arbitrary); `++svd-vals` for the values alone.  Needs
  `rows >= cols` — transpose a wide matrix and swap `u`/`v`.
- Vector helpers, on rank-1 rays of shape `~[n]` (as `++diag` returns, not
  `n x 1` matrices): `++matvec`, `++dotv`, `++nrm2`, `++axpyv`, `++col-dot`.

Every inner product in these arms sums LEFT TO RIGHT from the kind's `+0` — the
order a C or Rust kernel walks, so they stay jettable, and not NumPy's pairwise
summation.

Set rounding mode and tolerance with `++sake` before calling `++eig`, `++cg`,
`++pcg` or `++svd` (the bare `++sa` default `rtol` is unusable); `rtol`'s width
must match the component.

Random array filling (`rand-spec.md` section 8, in `/Users/neal/urbit/numerics/librand/`):

- `++fill-uniform`, `++fill-normal`, `++fill-expon` — fill an `%i754` ray
  (bloq 5/6, `@rs`/`@rd` only) element-by-element in row-major order from
  an `/lib/rand` engine. Philox elements get per-element counter
  treatment so a future jet can parallelize across elements and land
  identical bits regardless of thread scheduling: `++fill-uniform`'s
  single non-rejecting draw assigns element `i` counter `ctr0+i`
  directly; the rejection-based `++fill-normal`/`++fill-expon` instead
  give each element a `ctr0 + i*2^32` counter *window* to walk freely
  within, with an explicit crash if one element's rejection loop ever
  exhausts its own window.
- `++fill-below` — fill a `%uint` ray via Lemire's unbiased method, same
  per-element windowing.
- Posit (`%unum`) rays are deferred (see librand's `NEXT-STEPS.md`).

##  References

- Milton Abramowitz & Irene Stegun, _Handbook of Mathematical Functions with Formulas, Graphs, and Mathematical Tables_.  1964–2010.
- Forman Acton, _Numerical Methods that (Usually) Work_, 1ed.  1997.
- [Bartosz Ciechanowski, “Float Exposed” (webapp)](https://float.exposed/0x00000001)
- [David Goldberg, “What Every Computer Scientist Should Know About Floating-Point Arithmetic”](https://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html)
- Parviz Moin, _Fundamentals of Engineering Numerical Analysis_. 2ed.  2001.
