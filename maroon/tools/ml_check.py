#!/usr/bin/env python3
"""Oracle for /lib/maroon's PCA and k-means.

Checks the expectations asserted by

    maroon/desk/tests/lib/maroon-pca.hoon
    maroon/desk/tests/lib/maroon-kmeans.hoon

Run with no arguments; exits nonzero on any mismatch.

The Hoon semantics mirrored here, exactly:

  * a dataset is n x d, one row per sample;
  * +center subtracts the column means; +pca decomposes the CENTERED data with
    Saloon's one-sided Jacobi SVD (not an eigendecomposition of the
    covariance), keeps the top k right-singular vectors as the columns of
    .comp, and reports variances s_i^2/(n-1);
  * +assign is argmin over squared distances, ties to the lower centroid;
  * +update is the mean of each cluster, an empty cluster keeping its previous
    centroid;
  * +kmeans iterates to the exact fixed point of the assignment (not a
    tolerance), reporting the iteration count;
  * every sum folds LEFT TO RIGHT from +0.

Exact cases are verified op-by-op in exact rational arithmetic, so the Hoon can
compare whole rays; approximate ones come from numpy/sklearn.
"""

import sys
from fractions import Fraction

import numpy as np
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

FAIL = []
N = [0]


def check(label, got, want):
    N[0] += 1
    if got != want:
        FAIL.append(f"{label}: got {got!r}, want {want!r}")
        print(f"  FAIL {label}: got {got!r}, want {want!r}")


def close(label, got, want, tol=1e-9):
    N[0] += 1
    g, w = np.asarray(got, dtype=float), np.asarray(want, dtype=float)
    if g.shape != w.shape or not np.allclose(g, w, rtol=tol, atol=tol):
        FAIL.append(f"{label}: got {g.tolist()!r}, want {w.tolist()!r}")
        print(f"  FAIL {label}: got {g.tolist()!r}, want {w.tolist()!r}")


# ------------------------------------------------------- exact rational layer

def qmean(col):
    return sum((Fraction(v) for v in col), Fraction(0)) / len(col)


def qcenter(x):
    cols = list(zip(*x))
    means = [qmean(c) for c in cols]
    return [[Fraction(v) - m for v, m in zip(row, means)] for row in x], means


def exact_repr(fr):
    """A Fraction that a binary float represents exactly, as a short decimal."""
    f = float(fr)
    assert Fraction(f) == fr, f"{fr} is not exactly representable"
    return f


print("== +col-mean / +center: exact ==")
X1 = [[1, 0], [2, 0], [3, 0]]
xc, means = qcenter(X1)
check("col-mean [[1,0],[2,0],[3,0]]", [exact_repr(m) for m in means], [2.0, 0.0])
check("center   [[1,0],[2,0],[3,0]]",
      [[exact_repr(v) for v in r] for r in xc], [[-1.0, 0.0], [0.0, 0.0], [1.0, 0.0]])

X2 = [[0, 0], [4, 0], [0, 2], [4, 2]]
xc2, means2 = qcenter(X2)
check("col-mean 4x2", [exact_repr(m) for m in means2], [2.0, 1.0])
check("center 4x2", [[exact_repr(v) for v in r] for r in xc2],
      [[-2.0, -1.0], [2.0, -1.0], [-2.0, 1.0], [2.0, 1.0]])

print("== +pca: the axis is exact; the variance goes through sqrt-then-square ==")
#  X1 centered is [[-1,0],[0,0],[1,0]]: columns already orthogonal, so the
#  one-sided Jacobi does no rotation at all and V is the identity.
p = PCA(n_components=1).fit(np.array(X1, dtype=float))
close("pca X1 component (numpy)", np.abs(p.components_[0]), [1.0, 0.0])
close("pca X1 variance (numpy)", p.explained_variance_, [1.0])
#  what the Hoon actually returns for the variance: (sqrt(2))^2 / 2 in float64
s = np.sqrt(np.float64(2.0))
hoon_var = np.float64(np.float64(s * s) / np.float64(2.0))
print(f"  +pca variance for X1 = {hoon_var!r}  (1.0 to within an ulp)")
check("pca X1 variance is 1 ulp low", hoon_var == 1.0, False)
close("pca X1 variance ~ 1", [hoon_var], [1.0], tol=1e-15)

#  a case with a rotated axis, checked against sklearn
X3 = [[2.0, 1.0], [4.0, 3.0], [6.0, 5.0], [8.0, 7.0]]
p3 = PCA(n_components=2).fit(np.array(X3))
print(f"  pca X3 components (sklearn): {p3.components_.tolist()}")
print(f"  pca X3 variances  (sklearn): {p3.explained_variance_.tolist()}")
close("pca X3 first component is the 45-degree axis",
      np.abs(p3.components_[0]), [np.sqrt(0.5), np.sqrt(0.5)])
close("pca X3 variances", p3.explained_variance_, [20.0 / 3.0 * 2, 0.0], tol=1e-8)
#  scores of the training data on the first component, sign-free
scores = PCA(n_components=1).fit_transform(np.array(X3))
print(f"  pca X3 scores (sklearn): {np.abs(scores).ravel().tolist()}")

print("== +cov: exact at ddof=0 for this data ==")
#  Xc^T*Xc / n for X2, whose centred columns are orthogonal
cols = list(zip(*xc2))
gram = [[sum((a * b for a, b in zip(ci, cj)), Fraction(0)) for cj in cols] for ci in cols]
cov0 = [[v / len(X2) for v in r] for r in gram]
check("cov X2 ddof=0", [[exact_repr(v) for v in r] for r in cov0], [[4.0, 0.0], [0.0, 1.0]])
#  ddof=1 divides by 3 and is NOT exact, which is why the tests use ddof=0
N[0] += 1
if Fraction(float(gram[0][0] / 3)) == gram[0][0] / 3:
    FAIL.append("cov X2 ddof=1 unexpectedly exact; the test design assumes it is not")

print("== +assign / +update / +kmeans: exact ==")
#  1-D, four points, two obvious clusters
P = [[0.0], [1.0], [10.0], [11.0]]
C0 = [[0.0], [10.0]]


def assign(x, c):
    out = []
    for row in x:
        best, bd = 0, None
        for j, cen in enumerate(c):
            d = sum((Fraction(a) - Fraction(b)) ** 2 for a, b in zip(row, cen))
            if bd is None or d < bd:           # strict: ties keep the lower index
                best, bd = j, d
        out.append(best)
    return out


def update(x, lab, c, k):
    out = [list(map(Fraction, row)) for row in c]
    for j in range(k):
        rows = [row for row, l in zip(x, lab) if l == j]
        if not rows:
            continue
        out[j] = [sum((Fraction(r[t]) for r in rows), Fraction(0)) / len(rows)
                  for t in range(len(rows[0]))]
    return out


def kmeans(x, c0, maxit):
    c, lab, it = [list(map(Fraction, r)) for r in c0], assign(x, c0), 0
    while it < maxit:
        c2 = update(x, lab, c, len(c))
        lab2 = assign(x, c2)
        if lab2 == lab:
            return lab2, c2, it + 1
        c, lab, it = c2, lab2, it + 1
    return lab, c, it


check("assign P against C0", assign(P, C0), [0, 0, 1, 1])
check("update P", [[exact_repr(v) for v in r] for r in update(P, [0, 0, 1, 1], C0, 2)],
      [[0.5], [10.5]])
lab, c, it = kmeans(P, C0, 20)
check("kmeans P labels", lab, [0, 0, 1, 1])
check("kmeans P centroids", [[exact_repr(v) for v in r] for r in c], [[0.5], [10.5]])
check("kmeans P iterations", it, 1)
#  inertia: each point is 0.5 from its centroid, so 4 * 0.25
inertia = sum((Fraction(row[0]) - c[l][0]) ** 2 for row, l in zip(P, lab))
check("kmeans P inertia", exact_repr(inertia), 1.0)

#  an empty cluster keeps its previous centroid
check("update with an empty cluster",
      [[exact_repr(v) for v in r] for r in update(P, [0, 0, 0, 0], [[0.0], [99.0]], 2)],
      [[5.5], [99.0]])

#  2-D, two well-separated clusters, exact means
P2 = [[0.0, 0.0], [0.0, 2.0], [10.0, 0.0], [10.0, 2.0]]
lab2, c2, it2 = kmeans(P2, [[0.0, 0.0], [10.0, 0.0]], 20)
check("kmeans P2 labels", lab2, [0, 0, 1, 1])
check("kmeans P2 centroids", [[exact_repr(v) for v in r] for r in c2],
      [[0.0, 1.0], [10.0, 1.0]])

#  sklearn agrees on the partition (its own seeding, so compare the grouping)
km = KMeans(n_clusters=2, n_init=10, random_state=0).fit(np.array(P2))
groups = sorted(sorted(np.flatnonzero(km.labels_ == j).tolist()) for j in range(2))
check("kmeans P2 partition vs sklearn", groups, [[0, 1], [2, 3]])

print("== +linreg / +ridge / +predict / +mse / +r2: exact ==")
#  x = [0,0,2,2], y = 3x + 1.  Centred x is [-1,-1,1,1] with norm exactly 2, so
#  the one Householder step and the back-substitution stay exact.
RX = [[0], [0], [2], [2]]
RY = [1, 1, 7, 7]
xm = qmean([r[0] for r in RX])
ym = qmean(RY)
xc = [Fraction(r[0]) - xm for r in RX]
yc = [Fraction(v) - ym for v in RY]
#  OLS on one centred feature: coef = (xc.yc)/(xc.xc)
coef = sum((a * b for a, b in zip(xc, yc)), Fraction(0)) / sum((a * a for a in xc), Fraction(0))
check("linreg coef", exact_repr(coef), 3.0)
check("linreg intercept", exact_repr(ym - xm * coef), 1.0)
#  the Householder route the Hoon actually takes: norm 2, v = xc - alpha*e1
nx = Fraction(2)
alpha = nx if xc[0] < 0 else -nx
v = [xc[0] - alpha] + xc[1:]
vv = sum((t * t for t in v), Fraction(0))
f_q = (2 * sum((t * e for t, e in zip(v, [1, 0, 0, 0])), Fraction(0))) / vv
check("linreg Householder scale is exact", Fraction(float(f_q)) == f_q, True)
#  ridge with alpha = 12: (xc.xc + 12) coef = xc.yc -> 16 coef = 12
rc = sum((a * b for a, b in zip(xc, yc)), Fraction(0)) / (sum((a * a for a in xc), Fraction(0)) + 12)
check("ridge coef (alpha=12)", exact_repr(rc), 0.75)
check("ridge intercept (alpha=12)", exact_repr(ym - xm * rc), 3.25)
preds = [ym - xm * rc + rc * Fraction(r[0]) for r in RX]
res = [Fraction(yv) - p for yv, p in zip(RY, preds)]
mse = sum((t * t for t in res), Fraction(0)) / len(res)
check("ridge mse", exact_repr(mse), 5.0625)
sst = sum((t * t for t in yc), Fraction(0))
r2 = 1 - sum((t * t for t in res), Fraction(0)) / sst
check("ridge r2", exact_repr(r2), 0.4375)
#  and the OLS fit is perfect
check("linreg mse is 0", exact_repr(sum(((Fraction(yv) - (ym - xm * coef + coef * Fraction(r[0]))) ** 2
                                          for yv, r in zip(RY, RX)), Fraction(0))), 0.0)
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.metrics import mean_squared_error, r2_score
lr = LinearRegression().fit(np.array(RX, dtype=float), np.array(RY, dtype=float))
close("linreg vs sklearn", [lr.coef_[0], lr.intercept_], [3.0, 1.0])
rg = Ridge(alpha=12.0).fit(np.array(RX, dtype=float), np.array(RY, dtype=float))
close("ridge vs sklearn", [rg.coef_[0], rg.intercept_], [0.75, 3.25])
close("mse vs sklearn", [mean_squared_error(RY, rg.predict(np.array(RX, dtype=float)))], [5.0625])
close("r2 vs sklearn", [r2_score(RY, rg.predict(np.array(RX, dtype=float)))], [0.4375])

print("== +linreg / +ridge: two features, against sklearn ==")
MX = np.array([[1.0, 2.0], [2.0, 1.0], [3.0, 4.0], [4.0, 3.0], [5.0, 6.0]])
MY = np.array([3.1, 2.9, 7.2, 6.8, 11.1])
lr2 = LinearRegression().fit(MX, MY)
rg2 = Ridge(alpha=1.0).fit(MX, MY)
print(f"  linreg 2-feature: coef={lr2.coef_.tolist()} intercept={lr2.intercept_!r}")
print(f"  ridge  2-feature: coef={rg2.coef_.tolist()} intercept={rg2.intercept_!r}")
print(f"  linreg r2 on its training data: {r2_score(MY, lr2.predict(MX))!r}")

print("== ties go to the lower centroid index ==")
#  a point equidistant from both centroids
check("assign tie", assign([[5.0]], [[0.0], [10.0]]), [0])

print()
print(f"{N[0]} checks, {len(FAIL)} failures")
for f in FAIL:
    print(f"  {f}")
sys.exit(1 if FAIL else 0)
