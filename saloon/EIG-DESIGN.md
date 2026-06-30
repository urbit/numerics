# Design: Saloon eigendecomposition (`eig`)

Status: **A1 (symmetric Jacobi) and A2 (Hermitian Jacobi) shipped (PR #47); Phase B (Hessenberg + QR for general real matrices) pending.** 2026-06-06. The marquee consumer of Lagoon `%cplx`
(PR #46): eigenvalues/eigenvectors. Pure Hoon on the Lagoon `%i754`/`%cplx`
arrays; jets (SoftBLAS) only accelerate later, so this is fully parallel to the
SoftBLAS session.

Saloon (`saloon/desk/lib/saloon.hoon`) is a `sa` door over Lagoon's `la`,
carrying `rnd` (rounding) and `rtol` (relative tolerance). It already has the
scalar transcendentals (`sqrt`, `exp`, …) and pass-through matrix ops. `eig`
adds a numerical-linear-algebra layer on top.

---

## 1. Scope and the central split

A real square matrix has, in general, **complex** eigenvalues (in
conjugate pairs) — that is the whole reason `%cplx` had to exist first. But the
*symmetric/Hermitian* case is much easier and very common (covariance/Gram
matrices, graph Laplacians, quadratic forms), has **real** eigenvalues and
orthonormal eigenvectors, and needs no complex arithmetic at all.

So `eig` is phased:

- **Phase A — symmetric/Hermitian, via the Jacobi eigenvalue algorithm.**
  Real eigenvalues, orthonormal eigenvectors. Pure `%i754`. Robust, ~10 lines
  of rotation math, converges quadratically. A complete, useful capability on
  its own and the recommended first deliverable.
- **Phase B — general (nonsymmetric) real matrices, via Hessenberg reduction +
  Francis double-shift QR**, producing complex eigenvalues/eigenvectors on
  `%cplx`. Substantially harder; depends on `%cplx` and (for speed) its jets.

This doc specifies Phase A in full and sketches Phase B.

---

## 2. API

```
++  eigvals  |=(a=ray ...)   ::  -> ray  (1-D, the eigenvalues)
++  eigvecs  |=(a=ray ...)   ::  -> ray  (square, columns = eigenvectors)
++  eig      |=(a=ray ...)   ::  -> [vals=ray vecs=ray]
```

Phase A returns `%i754` reals (eigenvalues) and eigenvectors as columns:
`%i754` orthogonal matrix for symmetric (`%i754`) input; `%cplx` unitary
matrix for Hermitian (`%cplx`) input (A2 is shipped). Phase B returns `%cplx`
for general real.
A `?>` asserts squareness; Phase A also documents that it *assumes* symmetry
(it operates on the lower/upper triangle; a non-symmetric input is silently
symmetrized or rejected — see §4).

---

## 3. Building blocks to add to Saloon

Most are thin compositions of existing Lagoon arms (`transpose`, `mmul`, `dot`/
`dotc`, `sub`, `eye`, `get-item`, `set-item`, `scale`) and Saloon `sqrt`:

- `norm2` — Euclidean norm of a vector ray: `sqrt(Σ xᵢ²)` = `sqrt:sa` of
  `dot`(x,x) (real) / `dotc` (complex). 
- `normalize` — `x / norm2(x)`.
- `dagger` — conjugate transpose: `(conj (transpose a))` (= `transpose` for real).
- `identity` / `diag-of` / `off-diag-norm` — helpers for convergence tests.
- (Phase B) `householder` reflector and `givens` rotation builders.

These live in a new `+|  %linalg` section of `sa`.

---

## 4. Phase A — Jacobi (symmetric/Hermitian)

Classic cyclic Jacobi. Maintain `A` (working copy, converges to diagonal) and
`V` (accumulated rotations, converges to the eigenvector matrix, init `eye`).

Repeat sweeps until the off-diagonal Frobenius norm < `rtol·‖A‖` (or a sweep
cap 60 is hit — **emits a `~&` terminal trace (not a persistent log)**, never silently return a
non-converged result):

For each off-diagonal `(p,q)`, `p<q`, with `a_pq ≠ 0`, build a Givens rotation
that zeros `a_pq`. The angle needs only `sqrt` (no trig):

```
θ = (a_qq − a_pp) / (2·a_pq)
t = sign(θ) / (|θ| + sqrt(θ²+1))          ::  smaller root, for stability
c = 1 / sqrt(t²+1)
s = t·c
```

Apply the rotation on both sides: `A ← Jᵀ·A·J`, `V ← V·J`, where `J` is identity
except `J_pp=J_qq=c`, `J_pq=s`, `J_qp=−s`. Implement the row/column updates
directly with `get-item`/`set-item` on the four affected entries per `(p,q)`
(O(n) per rotation) rather than full `mmul` (O(n³)) — the scalar component
arithmetic uses the `rs`/`rd` door at the array's `bloq`.

Eigenvalues = `diag(A)` at convergence; eigenvectors = columns of `V`. (Eigenvalue order is unspecified (Jacobi pivot order); callers must sort
externally if needed.)

**Hermitian (`%cplx`) variant:** same skeleton with complex Givens rotations
and `dagger` instead of transpose; defer until after the real case works.

**Symmetry handling:** Phase A is only valid for symmetric input. Decision
(for review): (a) `?>` assert `is-close` symmetry and crash otherwise, or
(b) symmetrize `(A+Aᵀ)/2` silently, or (c) document "caller's responsibility."
Recommend **(a)** — assert, so a nonsymmetric matrix doesn't quietly get the
wrong algorithm (it should go to Phase B).

---

## 5. Phase B — general real → complex (sketch)

1. **Balancing** (optional, improves conditioning).
2. **Hessenberg reduction** via Householder reflectors: `A → H` upper-Hessenberg,
   orthogonally similar.
3. **Francis double-shift QR** iteration on `H`: implicit double-shift to keep
   real arithmetic until 1×1 (real eigenvalue) or 2×2 (complex-conjugate pair)
   blocks deflate off the diagonal. Complex eigenvalues come from the 2×2
   blocks' characteristic roots → `%cplx`.
4. **Eigenvectors** by back-substitution on the (quasi-)triangular factor, then
   transform back through the accumulated orthogonal/Householder transforms.

Depends on `%cplx` arithmetic (PR #46) and benefits from the `cgemm`/`zgemm`
jets (SoftBLAS session). This is the large follow-on; design it in full when
Phase A lands.

---

## 6. Numerical notes

- Saloon's transcendentals are naive, but `sqrt` (the only special function
  Jacobi needs) is fine; component arithmetic is exact IEEE via the float door.
- Convergence: iterate to `rtol` (Saloon already carries it); **emit a `~&`
  terminal trace (not a persistent log) if the sweep/iteration cap is hit**
  rather than returning a silent non-converged answer.
- Determinism: Jacobi is deterministic; results are reproducible.

---

## 7. Oracle & tests

NumPy: `numpy.linalg.eigh` (symmetric, Phase A) and `eig` (general, Phase B).
Caveats for cross-checking: eigenvectors are unique only up to **sign/phase**
and **ordering**, and degenerate eigenvalues give an arbitrary basis of the
eigenspace — so tests should check **invariants**, not raw vector equality:

- `A·v = λ·v` (residual `‖A·v − λ·v‖ < tol`) for each pair.
- `V` orthonormal: `Vᵀ·V ≈ I`.
- reconstruction `V·diag(λ)·Vᵀ ≈ A`.
- eigenvalue *set* matches `eigh` (sorted), within tol.

Curated small cases (2×2, 3×3 symmetric with known spectra, e.g.
`[[2,1],[1,2]] → {1,3}`), plus a `tools/eig_check.py` NumPy oracle in the
`saloon` tools dir.

---

## 8. Phasing

- **A1 (DONE, PR #47):** real symmetric Jacobi `eigvals`/`eigvecs`/`eig`, with
  invariant-based tests and the NumPy `eigh` oracle.
- **A2 (DONE):** Hermitian (`%cplx`) Jacobi. `eig` now dispatches on kind
  (`%cplx` → Hermitian, else symmetric). The complex rotation `J` that zeros
  `a_pq` has real diagonal `c` and complex off-diagonal `b = s·(a_pq/|a_pq|)`
  with `J[q,p] = −conj(b)`; updates are `A ← Jᴴ·A·J`, `V ← V·J`. Returns real
  `%i754` eigenvalues and `%cplx` (unitary) eigenvectors, matching NumPy
  `eigh`. Tested over `@cs`/`@cd` against `eigvalsh` (`tools/eigh_check.py`).
  Two robustness notes from this work: (1) the internal Newton `sqrt` is now
  iteration-capped (`fsqt`) — a sub-ULP `rtol` previously made it oscillate
  forever and once OOM-crashed the ship; (2) `rtol`'s width must match the
  component width (`@rs` for `@cs`, `@rd` for `@cd`); and (3) The Hermitian
  check uses `+cnear` (magnitude-based), so ±0.0 sign differences in conjugate
  pairs and rounding-induced mismatches within `+stol` are accepted.
- **B:** general real → complex via Hessenberg + double-shift QR (own design +
  PR), once `%cplx` (PR #46) has landed and ideally its jets exist.
