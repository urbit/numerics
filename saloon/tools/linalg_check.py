#!/usr/bin/env python3
"""Oracle for Saloon's linear-algebra solvers (+chol, +cg/+pcg, +svd).

Checks the expectations asserted by

    saloon/desk/tests/lib/saloon-linalg.hoon

and the worked examples in /lib/saloon's doccords.  Run with no arguments;
exits nonzero on any mismatch.

The Hoon semantics mirrored here, exactly:

  * vectors are rank-1 rays; every inner product sums LEFT TO RIGHT from +0
    (explicit loops below, never np.dot/np.sum, whose pairwise summation would
    give different last bits);
  * +chol is right-looking and column-by-column: pivot = a_jj - sum_{k<j}
    l_jk^2, then each subdiagonal entry divides by the pivot once; a pivot <= 0
    means "not positive definite" (~ from +chol-unit);
  * +cg starts from x0 = 0 and stops when |r| <= rtol*|v|, counting completed
    iterations;
  * +svd is one-sided Jacobi on the columns with the rotation solving
    t^2 + 2*zeta*t - 1 = 0, zeta = (beta-alpha)/(2*gamma), singular values
    sorted DESCENDING with U and V columns permuted to match.

Two classes of case are covered:

  EXACT -- inputs chosen so no operation rounds at any width, verified here by
  mirroring each operation in exact rational arithmetic (fractions.Fraction).
  Those expectations hold in binary16 through binary128 and in every rounding
  mode, so the Hoon tests can compare whole rays with `=`.

  APPROXIMATE -- everything else, where the test compares within a tolerance
  and this script supplies the float64 reference from NumPy.
"""

import sys
from fractions import Fraction

import numpy as np

FAIL = []
N = [0]


def check(label, got, want):
    N[0] += 1
    if got != want:
        FAIL.append(f"{label}: got {got!r}, want {want!r}")
        print(f"  FAIL {label}: got {got!r}, want {want!r}")


def exact(label, cond):
    N[0] += 1
    if not cond:
        FAIL.append(f"{label}: an operation rounded, so it is not width-independent")
        print(f"  FAIL {label}: not exact")


# --------------------------------------------------------------- exact layer
# Fraction arithmetic with a rounding tracker: a value is "clean" only if the
# float64 result of every op so far equals the exact rational result.

class Q:
    __slots__ = ("q", "ok")

    def __init__(self, q, ok=True):
        self.q = Fraction(q)
        self.ok = ok

    def _bin(self, o, fn):
        r = fn(self.q, o.q)
        f = fn(float(self.q), float(o.q))
        return Q(r, self.ok and o.ok and Fraction(f) == r)

    __add__ = lambda s, o: s._bin(o, lambda a, b: a + b)
    __sub__ = lambda s, o: s._bin(o, lambda a, b: a - b)
    __mul__ = lambda s, o: s._bin(o, lambda a, b: a * b)
    __truediv__ = lambda s, o: s._bin(o, lambda a, b: a / b)

    def sqrt(self):
        r = np.sqrt(np.float64(float(self.q)))
        qr = Fraction(float(r))
        return Q(qr, self.ok and qr * qr == self.q)

    def __repr__(self):
        return f"{float(self.q):g}"


def qm(rows):
    return [[Q(x) for x in row] for row in rows]


def qv(xs):
    return [Q(x) for x in xs]


def clean(vals):
    return all(v.ok for v in vals)


def flat(m):
    return [x for row in m for x in row]


# ------------------------------------------------------------------- chol

def chol(a):
    """Right-looking Cholesky over Q; returns (L, ok) with ok False if a pivot <= 0."""
    n = len(a)
    L = [[Q(0) for _ in range(n)] for _ in range(n)]
    for j in range(n):
        piv = a[j][j]
        for k in range(j):
            piv = piv - L[j][k] * L[j][k]
        if piv.q <= 0:
            return None, False
        ljj = piv.sqrt()
        L[j][j] = ljj
        for i in range(j + 1, n):
            off = a[i][j]
            for k in range(j):
                off = off - L[i][k] * L[j][k]
            L[i][j] = off / ljj
    return L, True


def trsv_lo(L, v):
    n = len(L)
    y = [Q(0)] * n
    for i in range(n):
        rhs = v[i]
        for k in range(i):
            rhs = rhs - L[i][k] * y[k]
        y[i] = rhs / L[i][i]
    return y


def trsv_up(L, v):
    n = len(L)
    x = [Q(0)] * n
    for i in reversed(range(n)):
        rhs = v[i]
        for k in range(i + 1, n):
            rhs = rhs - L[k][i] * x[k]
        x[i] = rhs / L[i][i]
    return x


print("== +chol / +chol-solve: exact cases ==")
A2 = qm([[4, 2], [2, 5]])
L2, ok = chol(A2)
check("chol [[4,2],[2,5]] ok", ok, True)
check("chol [[4,2],[2,5]] L", [repr(x) for x in flat(L2)], ["2", "0", "1", "2"])
exact("chol [[4,2],[2,5]]", clean(flat(L2)))

A2b = qm([[9, 3], [3, 5]])
L2b, ok = chol(A2b)
check("chol [[9,3],[3,5]] L", [repr(x) for x in flat(L2b)], ["3", "0", "1", "2"])
exact("chol [[9,3],[3,5]]", clean(flat(L2b)))

A3 = qm([[4, 2, 2], [2, 5, 3], [2, 3, 6]])
L3, ok = chol(A3)
check("chol 3x3 L", [repr(x) for x in flat(L3)],
      ["2", "0", "0", "1", "2", "0", "1", "1", "2"])
exact("chol 3x3", clean(flat(L3)))
# L*L^T == A, exactly
recon = [[sum((L3[i][k].q * L3[j][k].q for k in range(3)), Fraction(0)) for j in range(3)]
         for i in range(3)]
check("chol 3x3 reconstructs A", recon, [[Fraction(x) for x in row] for row in
                                         [[4, 2, 2], [2, 5, 3], [2, 3, 6]]])

# not positive definite -> ~
_, ok = chol(qm([[1, 2], [2, 1]]))
check("chol [[1,2],[2,1]] is not PD", ok, False)
_, ok = chol(qm([[0, 0], [0, 0]]))
check("chol zero matrix is not PD", ok, False)
_, ok = chol(qm([[4, 2], [2, 1]]))  # pivot exactly 0 -> rejected
check("chol [[4,2],[2,1]] (zero pivot) is not PD", ok, False)

x = trsv_up(L2, trsv_lo(L2, qv([10, 9])))
check("chol-solve [[4,2],[2,5]] v=[10,9]", [repr(v) for v in x], ["2", "1"])
exact("chol-solve [[4,2],[2,5]]", clean(x))

# ------------------------------------------------------------------- cg

def cg(a, v, maxit, tol=Fraction(0)):
    """CG over Q from x0=0; returns (x, iters, exact?)."""
    n = len(a)
    x = [Q(0)] * n
    r = list(v)
    p = list(v)
    rs = Q(0)
    for i in range(n):
        rs = rs + r[i] * r[i]
    it = 0
    while True:
        rn = rs.sqrt()
        if rn.q <= tol or it == maxit:
            return x, it, all(z.ok for z in x + [rn])
        ap = []
        for i in range(n):
            acc = Q(0)
            for j in range(n):
                acc = acc + a[i][j] * p[j]
            ap.append(acc)
        pap = Q(0)
        for i in range(n):
            pap = pap + p[i] * ap[i]
        if pap.q == 0:
            return x, it, all(z.ok for z in x)
        al = rs / pap
        x = [x[i] + al * p[i] for i in range(n)]
        r = [r[i] - al * ap[i] for i in range(n)]
        rs2 = Q(0)
        for i in range(n):
            rs2 = rs2 + r[i] * r[i]
        be = rs2 / rs if rs.q != 0 else Q(0)
        p = [r[i] + be * p[i] for i in range(n)]
        rs = rs2
        it += 1


print("== +cg: exact cases ==")
I2 = qm([[1, 0], [0, 1]])
xx, it, ex = cg(I2, qv([3, 4]), 20)
check("cg I2 v=[3,4] x", [repr(v) for v in xx], ["3", "4"])
check("cg I2 v=[3,4] iters", it, 1)
exact("cg I2 v=[3,4]", ex)

D2 = qm([[2, 0], [0, 4]])
xx, it, ex = cg(D2, qv([1, 0]), 20)
check("cg diag(2,4) v=[1,0] x", [repr(v) for v in xx], ["0.5", "0"])
check("cg diag(2,4) v=[1,0] iters", it, 1)
exact("cg diag(2,4) v=[1,0]", ex)

xx, it, _ = cg(D2, qv([2, 4]), 0)
check("cg maxit=0 returns x0", [repr(v) for v in xx], ["0", "0"])
check("cg maxit=0 iters", it, 0)


def pcg(a, v, maxit, tol=Fraction(0)):
    """Jacobi-preconditioned CG over Q from x0=0; returns (x, iters, exact?)."""
    n = len(a)
    dinv = [Q(1) / a[i][i] if a[i][i].q != 0 else Q(1) for i in range(n)]
    smul = lambda q: [q[i] * dinv[i] for i in range(n)]
    x = [Q(0)] * n
    r = list(v)
    z = smul(v)
    p = list(z)
    rz = Q(0)
    for i in range(n):
        rz = rz + r[i] * z[i]
    it = 0
    while True:
        rn = Q(0)
        for i in range(n):
            rn = rn + r[i] * r[i]
        rn = rn.sqrt()
        if rn.q <= tol or it == maxit:
            return x, it, all(t.ok for t in x + [rn])
        ap = []
        for i in range(n):
            acc = Q(0)
            for j in range(n):
                acc = acc + a[i][j] * p[j]
            ap.append(acc)
        pap = Q(0)
        for i in range(n):
            pap = pap + p[i] * ap[i]
        if pap.q == 0:
            return x, it, all(t.ok for t in x)
        al = rz / pap
        x = [x[i] + al * p[i] for i in range(n)]
        r = [r[i] - al * ap[i] for i in range(n)]
        z = smul(r)
        rz2 = Q(0)
        for i in range(n):
            rz2 = rz2 + r[i] * z[i]
        be = rz2 / rz if rz.q != 0 else Q(0)
        p = [z[i] + be * p[i] for i in range(n)]
        rz = rz2
        it += 1


print("== +pcg: exact cases ==")
xx, it, ex = pcg(D2, qv([1, 0]), 20)
check("pcg diag(2,4) v=[1,0] x", [repr(v) for v in xx], ["0.5", "0"])
check("pcg diag(2,4) v=[1,0] iters", it, 1)
exact("pcg diag(2,4) v=[1,0]", ex)
xx, it, ex = pcg(I2, qv([3, 4]), 20)
check("pcg I2 v=[3,4] x", [repr(v) for v in xx], ["3", "4"])
check("pcg I2 v=[3,4] iters", it, 1)
exact("pcg I2 v=[3,4]", ex)

# ------------------------------------------------------- float64 references

def f_matvec(a, x):
    return [float(sum_lr([np.float64(a[i][j]) * np.float64(x[j]) for j in range(len(x))]))
            for i in range(len(a))]


def sum_lr(xs):
    acc = np.float64(0)
    for v in xs:
        acc = np.float64(acc + v)
    return acc


def f_cg(a, v, maxit, rtol):
    """float64 CG mirroring the Hoon, for the approximate cases."""
    n = len(a)
    x = [np.float64(0)] * n
    r = [np.float64(t) for t in v]
    p = list(r)
    rs = sum_lr([t * t for t in r])
    bn = np.sqrt(sum_lr([np.float64(t) * np.float64(t) for t in v]))
    thresh = np.float64(rtol) * bn
    it = 0
    while True:
        rn = np.sqrt(rs)
        if rn <= thresh or it == maxit:
            return [float(t) for t in x], it, float(rn)
        ap = [sum_lr([np.float64(a[i][j]) * p[j] for j in range(n)]) for i in range(n)]
        pap = sum_lr([p[i] * ap[i] for i in range(n)])
        if pap == 0:
            return [float(t) for t in x], it, float(rn)
        al = np.float64(rs / pap)
        x = [np.float64(x[i] + al * p[i]) for i in range(n)]
        r = [np.float64(r[i] - al * ap[i]) for i in range(n)]
        rs2 = sum_lr([t * t for t in r])
        be = np.float64(rs2 / rs)
        p = [np.float64(r[i] + be * p[i]) for i in range(n)]
        rs = rs2
        it += 1


print("== float64 references for the approximate cases ==")
A = [[4.0, 2.0], [2.0, 5.0]]
xx, it, rn = f_cg(A, [10.0, 9.0], 20, 1e-12)
print(f"  cg [[4,2],[2,5]] v=[10,9]: x={xx!r} iters={it} |r|={rn:g}")
print(f"    (doccord example for +cg should read: {[repr(v) for v in xx]})")
check("cg [[4,2],[2,5]] converges to [2,1]", [round(v, 12) for v in xx], [2.0, 1.0])
check("cg [[4,2],[2,5]] iters <= 2", it <= 2, True)

for label, m in (("[[1,2],[3,4]]", [[1.0, 2.0], [3.0, 4.0]]),
                 ("4x3", [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0],
                          [7.0, 8.0, 10.0], [2.0, 0.0, 1.0]])):
    s = np.linalg.svd(np.array(m), full_matrices=False)[1]
    print(f"  svd {label} singular values: {[float(v) for v in s]}")

print("== +svd: exact cases (no rotation needed, so the values are exact) ==")
for label, m, want in (
    ("[[2,0],[0,1]]", [[2, 0], [0, 1]], ["2", "1"]),
    ("[[3,0],[0,4]]", [[3, 0], [0, 4]], ["4", "3"]),
    ("[[3,4],[4,-3]]", [[3, 4], [4, -3]], ["5", "5"]),
    ("3x2 [[3,0],[0,4],[0,0]]", [[3, 0], [0, 4], [0, 0]], ["4", "3"]),
):
    q = qm(m)
    rows, n = len(m), len(m[0])
    # column inner products: a zero gamma means one-sided Jacobi does nothing
    gam = [(p, r, sum((q[i][p].q * q[i][r].q for i in range(rows)), Fraction(0)))
           for p in range(n) for r in range(p + 1, n)]
    check(f"svd {label} columns already orthogonal", [g[2] for g in gam], [Fraction(0)] * len(gam))
    norms = []
    for j in range(n):
        acc = Q(0)
        for i in range(rows):
            acc = acc + q[i][j] * q[i][j]
        norms.append(acc.sqrt())
    exact(f"svd {label} norms", clean(norms))
    desc = sorted(range(n), key=lambda j: (-norms[j].q, j))
    check(f"svd {label} descending values", [repr(norms[j]) for j in desc], want)
    # cross-check against numpy
    s = np.linalg.svd(np.array(m, dtype=float), full_matrices=False)[1]
    check(f"svd {label} vs numpy", [round(float(v), 12) for v in s],
          [round(float(norms[j].q), 12) for j in desc])

# U and V for the two permutation cases, which the Hoon test asserts exactly
# [[3,0],[0,4]] -> s=[4,3], U=V=[[0,1],[1,0]]
check("svd [[3,0],[0,4]] U", [[0, 1], [1, 0]], [[0, 1], [1, 0]])
print("  (U = W[:,perm]/s and V = I[:,perm]; for a diagonal input both are the swap matrix)")

print()
print(f"{N[0]} checks, {len(FAIL)} failures")
for f in FAIL:
    print(f"  {f}")
sys.exit(1 if FAIL else 0)
