#!/usr/bin/env python3
"""Oracle for the axis/sort/distance arms of /lib/lagoon.

Recomputes every expectation asserted by

    lagoon/desk/tests/lib/lagoon-axis.hoon
    lagoon/desk/tests/lib/lagoon-sort-dist.hoon
    lagoon/desk/tests/lib/lagoon-axis-rounding.hoon

and checks it against the literal written there, so the test suite has a
reference outside Hoon.  Run with no arguments; exits nonzero on any mismatch.

The Hoon semantics mirrored here, exactly:

  * row-major storage; reducing dim d removes it from the shape (rank-1 input
    reduces to shape [1]);
  * every reduction folds LEFT TO RIGHT along the axis, seeded with the slice's
    first element -- NOT numpy's pairwise summation, so the folds below are
    explicit Python loops over numpy scalars of the target dtype;
  * mean = one division of the sequential sum; var = two-pass with one final
    division; std = one final sqrt;
  * sort/argsort break ties toward the lower original index (stable);
  * cdist-sq sums squared differences DIRECTLY, never via the Gram identity.

Widths: binary16/32/64 come from numpy.  binary128 is deliberately absent --
numpy's longdouble is not IEEE binary128 (on arm64 macOS it is plain float64),
so it would silently answer as f64.  Every test case that covers bloq 7 is
therefore built from values whose every intermediate is exact, which this
script verifies with exact rational arithmetic (`fractions.Fraction`): if no
operation rounds at any width, the result is the same in all four widths and in
all four rounding modes.
"""

import struct
import sys
from fractions import Fraction

import numpy as np

FAILURES = []
CHECKS = [0]


def check(label, got, want):
    CHECKS[0] += 1
    if got != want:
        FAILURES.append(f"{label}: got {got!r}, want {want!r}")
        print(f"  FAIL {label}: got {got!r}, want {want!r}")


# ---------------------------------------------------------------- exactness

def exact_ops_only(fn):
    """Run fn against a checking scalar type; return True if nothing rounded.

    fn is called with an arithmetic namespace and must return a list of
    numbers.  Each op is done in both float64 and exact rationals; a rounding
    anywhere (including the final sqrt) marks the case inexact.
    """
    rounded = [False]

    class Chk:
        __slots__ = ("f", "q")

        def __init__(self, f, q):
            self.f, self.q = f, q

        def _wrap(self, f, q):
            if Fraction(f) != q:
                rounded[0] = True
            return Chk(f, q)

        def __add__(self, o):
            return self._wrap(self.f + o.f, self.q + o.q)

        def __sub__(self, o):
            return self._wrap(self.f - o.f, self.q - o.q)

        def __mul__(self, o):
            return self._wrap(self.f * o.f, self.q * o.q)

        def __truediv__(self, o):
            return self._wrap(self.f / o.f, self.q / o.q)

        def sqrt(self):
            r = np.sqrt(np.float64(self.f))
            q = Fraction(float(r))
            if q * q != self.q:
                rounded[0] = True
            return Chk(float(r), q)

        def __lt__(self, o):
            return self.f < o.f

        def __gt__(self, o):
            return self.f > o.f

        def __eq__(self, o):
            return self.f == o.f

    vals = fn(lambda x: Chk(float(x), Fraction(x)))
    return (not rounded[0]), [v.f for v in vals]


# ------------------------------------------------------------ float helpers

DT = {16: np.float16, 32: np.float32, 64: np.float64}
FMT = {16: ("<e", "<H", 4), 32: ("<f", "<I", 8), 64: ("<d", "<Q", 16)}


def bits(width, x):
    packfmt, intfmt, digits = FMT[width]
    raw = struct.unpack(intfmt, struct.pack(packfmt, DT[width](x)))[0]
    return f"0x{raw:0{digits}x}"


def fold(dt, xs, op):
    """Left-to-right fold seeded with the head, in dtype dt."""
    acc = dt(xs[0])
    for x in xs[1:]:
        acc = dt(op(acc, dt(x)))
    return acc


def add_f(dt, xs):
    return fold(dt, xs, lambda a, b: a + b)


def slices(a, dim):
    """Every 1-D slice along dim, in the row-major order of the reduced shape."""
    a = np.asarray(a)
    moved = np.moveaxis(a, dim, -1)
    return moved.reshape(-1, a.shape[dim]).tolist(), list(moved.shape[:-1])


def reduce_dim(a, dim, dt, fn):
    sl, shape = slices(a, dim)
    out = [fn(dt, s) for s in sl]
    return [float(x) for x in out], (shape or [1])


def flat(a):
    return [float(x) for x in np.asarray(a).reshape(-1)]


# --------------------------------------------------- directed rounding (f32/f64)

PREC = {32: (24, -126), 64: (53, -1022)}


def round_mode(frac, width, mode):
    """Round an exact Fraction to binary<width> under mode n/z/u/d."""
    p, emin = PREC[width]
    if frac == 0:
        return Fraction(0)
    sign = -1 if frac < 0 else 1
    f = abs(frac)
    e = 0
    while Fraction(2) ** (e + 1) <= f:
        e += 1
    while Fraction(2) ** e > f:
        e -= 1
    q = max(e - (p - 1), emin - (p - 1))
    scaled = f / Fraction(2) ** q
    lo = scaled.numerator // scaled.denominator
    rem = scaled - lo
    if rem == 0:
        n = lo
    elif mode == "n":
        n = lo if rem < Fraction(1, 2) else lo + 1 if rem > Fraction(1, 2) else lo + (lo & 1)
    elif mode == "z":
        n = lo
    elif mode == "u":
        n = lo + 1 if sign > 0 else lo
    elif mode == "d":
        n = lo if sign > 0 else lo + 1
    else:
        raise ValueError(mode)
    return sign * n * Fraction(2) ** q


def sum_mode(xs, width, mode):
    """Left-to-right sum of exact Fractions, rounding each add under mode."""
    acc = round_mode(Fraction(xs[0]), width, mode)
    for x in xs[1:]:
        acc = round_mode(acc + round_mode(Fraction(x), width, mode), width, mode)
    return acc


def frac_bits(frac, width):
    return bits(width, float(frac))


# ============================================================== TASK: /lagoon-axis

A = [[1, 2, 3], [3, 4, 5]]
B = [[1, 3], [2, 6]]
C = [1, 3]
D = [3, 4]
E = [[3, 4], [6, 8]]
F = [[1, 2], [2, 4]]
CUB = [[[0, 1], [2, 3]], [[4, 5], [6, 7]]]

print("== lagoon-axis: reductions (expect exact in f16/f32/f64) ==")
for w, dt in DT.items():
    tag = f"f{w}"
    check(f"{tag} sum-dim A 0", reduce_dim(A, 0, dt, add_f)[0], [4, 6, 8])
    check(f"{tag} sum-dim A 1", reduce_dim(A, 1, dt, add_f)[0], [6, 12])
    check(f"{tag} sum-dim C 0", reduce_dim(C, 0, dt, add_f)[0], [4])
    check(f"{tag} sum-dim C shape", reduce_dim(C, 0, dt, add_f)[1], [1])
    mul_f = lambda d, xs: fold(d, xs, lambda a, b: a * b)
    check(f"{tag} prod-dim A 0", reduce_dim(A, 0, dt, mul_f)[0], [3, 8, 15])
    check(f"{tag} prod-dim A 1", reduce_dim(A, 1, dt, mul_f)[0], [6, 60])
    mx = lambda d, xs: fold(d, xs, lambda a, b: a if a > b else b)
    mn = lambda d, xs: fold(d, xs, lambda a, b: a if a < b else b)
    check(f"{tag} max-dim A 0", reduce_dim(A, 0, dt, mx)[0], [3, 4, 5])
    check(f"{tag} max-dim A 1", reduce_dim(A, 1, dt, mx)[0], [3, 5])
    check(f"{tag} min-dim A 0", reduce_dim(A, 0, dt, mn)[0], [1, 2, 3])
    check(f"{tag} min-dim A 1", reduce_dim(A, 1, dt, mn)[0], [1, 3])
    # arg*: first extremum
    amax = lambda d, xs: float(int(np.argmax(np.array(xs, dtype=d))))
    amin = lambda d, xs: float(int(np.argmin(np.array(xs, dtype=d))))
    check(f"{tag} argmax-dim A 0", reduce_dim(A, 0, dt, amax)[0], [1, 1, 1])
    check(f"{tag} argmax-dim A 1", reduce_dim(A, 1, dt, amax)[0], [2, 2])
    check(f"{tag} argmin-dim A 0", reduce_dim(A, 0, dt, amin)[0], [0, 0, 0])
    check(f"{tag} argmin-dim A 1", reduce_dim(A, 1, dt, amin)[0], [0, 0])
    check(f"{tag} argmax tie [5,5,1]", reduce_dim([5, 5, 1], 0, dt, amax)[0], [0])
    check(f"{tag} argmin tie [1,1,5]", reduce_dim([1, 1, 5], 0, dt, amin)[0], [0])
    # mean / var / std
    mean = lambda d, xs: d(add_f(d, xs) / d(len(xs)))
    check(f"{tag} mean-dim A 0", reduce_dim(A, 0, dt, mean)[0], [2, 3, 4])
    check(f"{tag} mean-dim A 1", reduce_dim(A, 1, dt, mean)[0], [2, 4])
    check(f"{tag} mean A whole", [float(mean(dt, flat(A)))], [3])
    check(f"{tag} mean C whole", [float(mean(dt, flat(C)))], [2])

    def var(d, xs, ddof=0):
        m = d(add_f(d, xs) / d(len(xs)))
        sq = [d(d(d(x) - m) * d(d(x) - m)) for x in xs]
        return d(add_f(d, sq) / d(len(xs) - ddof))

    check(f"{tag} var-dim B 1 ddof0", reduce_dim(B, 1, dt, lambda d, s: var(d, s, 0))[0], [1, 4])
    check(f"{tag} var-dim B 1 ddof1", reduce_dim(B, 1, dt, lambda d, s: var(d, s, 1))[0], [2, 8])
    check(f"{tag} std-dim B 1 ddof0",
          reduce_dim(B, 1, dt, lambda d, s: d(np.sqrt(var(d, s, 0))))[0], [1, 2])
    check(f"{tag} var C ddof0", [float(var(dt, flat(C), 0))], [1])
    check(f"{tag} var C ddof1", [float(var(dt, flat(C), 1))], [2])
    check(f"{tag} std C ddof0", [float(dt(np.sqrt(var(dt, flat(C), 0))))], [1])
    # norms
    l1 = lambda d, xs: add_f(d, [abs(d(x)) for x in xs])
    l2 = lambda d, xs: d(np.sqrt(add_f(d, [d(d(x) * d(x)) for x in xs])))
    li = lambda d, xs: fold(d, [abs(d(x)) for x in xs], lambda a, b: a if a > b else b)
    check(f"{tag} norm D l1", [float(l1(dt, D))], [7])
    check(f"{tag} norm D l2", [float(l2(dt, D))], [5])
    check(f"{tag} norm D linf", [float(li(dt, D))], [4])
    check(f"{tag} norm F fro", [float(l2(dt, flat(F)))], [5])
    check(f"{tag} norm [-3,4] l2", [float(l2(dt, [-3, 4]))], [5])
    check(f"{tag} norm-dim E 1 l2", reduce_dim(E, 1, dt, l2)[0], [5, 10])
    check(f"{tag} norm-dim E 1 l1", reduce_dim(E, 1, dt, l1)[0], [7, 14])
    check(f"{tag} norm-dim E 1 linf", reduce_dim(E, 1, dt, li)[0], [4, 8])
    check(f"{tag} norm-dim E 0 l1", reduce_dim(E, 0, dt, l1)[0], [9, 12])
    # rank 3
    check(f"{tag} sum-dim CUB 1", reduce_dim(CUB, 1, dt, add_f)[0], [2, 4, 10, 12])
    check(f"{tag} sum-dim CUB 0", reduce_dim(CUB, 0, dt, add_f)[0], [4, 6, 8, 10])
    check(f"{tag} sum-dim CUB 2", reduce_dim(CUB, 2, dt, add_f)[0], [1, 5, 9, 13])
    check(f"{tag} max-dim CUB 1", reduce_dim(CUB, 1, dt, mx)[0], [2, 3, 6, 7])
    check(f"{tag} min-dim CUB 1", reduce_dim(CUB, 1, dt, mn)[0], [0, 1, 4, 5])
    check(f"{tag} argmax-dim CUB 1", reduce_dim(CUB, 1, dt, amax)[0], [1, 1, 1, 1])
    check(f"{tag} mean-dim CUB 1", reduce_dim(CUB, 1, dt, mean)[0], [1, 2, 5, 6])

print("== lagoon-axis: broadcasting (numpy is the reference) ==")
check("broadcast row", flat(np.broadcast_to(np.array([[1, 2, 3]]), (2, 3))), [1, 2, 3, 1, 2, 3])
check("broadcast col", flat(np.broadcast_to(np.array([[1], [2]]), (2, 3))), [1, 1, 1, 2, 2, 2])
check("broadcast rank-up", flat(np.broadcast_to(np.array([1, 2, 3]), (2, 3))), [1, 2, 3, 1, 2, 3])
check("broadcast scalar", flat(np.broadcast_to(np.array([[7]]), (2, 2))), [7, 7, 7, 7])
check("broadcast identity", flat(np.broadcast_to(np.array(A), (2, 3))), flat(A))
check("broadcast rank-3 mid",
      flat(np.broadcast_to(np.array(CUB)[:, :1, :], (2, 3, 2))),
      [0, 1, 0, 1, 0, 1, 4, 5, 4, 5, 4, 5])
check("broadcast then add",
      flat(np.array(A) + np.broadcast_to(np.array([[1, 2, 3]]), (2, 3))),
      [2, 4, 6, 4, 6, 8])

print("== exactness: no operation may round, at any width ==")


def exact_case(label, fn):
    ok, vals = exact_ops_only(fn)
    CHECKS[0] += 1
    if not ok:
        FAILURES.append(f"{label}: an operation rounded, so it is NOT width-independent")
        print(f"  FAIL {label}: rounds")
    return vals


exact_case("sum A dim0", lambda mk: [mk(1) + mk(3), mk(2) + mk(4), mk(3) + mk(5)])
exact_case("prod A dim1", lambda mk: [mk(1) * mk(2) * mk(3), mk(3) * mk(4) * mk(5)])
exact_case("mean A dim0", lambda mk: [(mk(1) + mk(3)) / mk(2)])
exact_case("mean A whole", lambda mk: [(mk(1) + mk(2) + mk(3) + mk(3) + mk(4) + mk(5)) / mk(6)])
exact_case("var B row0 ddof0",
           lambda mk: [((mk(1) - mk(2)) * (mk(1) - mk(2)) + (mk(3) - mk(2)) * (mk(3) - mk(2))) / mk(2)])
exact_case("std B row1 ddof0",
           lambda mk: [(((mk(2) - mk(4)) * (mk(2) - mk(4)) + (mk(6) - mk(4)) * (mk(6) - mk(4))) / mk(2)).sqrt()])
exact_case("norm D l2", lambda mk: [(mk(3) * mk(3) + mk(4) * mk(4)).sqrt()])
exact_case("norm F fro",
           lambda mk: [(mk(1) * mk(1) + mk(2) * mk(2) + mk(2) * mk(2) + mk(4) * mk(4)).sqrt()])
exact_case("norm-dim E l2", lambda mk: [(mk(3) * mk(3) + mk(4) * mk(4)).sqrt(),
                                        (mk(6) * mk(6) + mk(8) * mk(8)).sqrt()])
exact_case("cdist pts", lambda mk: [(mk(1) - mk(0)) * (mk(1) - mk(0)) + (mk(0) - mk(1)) * (mk(0) - mk(1))])
exact_case("cdist 3-4-5", lambda mk: [(mk(0) - mk(3)) * (mk(0) - mk(3)) + (mk(0) - mk(4)) * (mk(0) - mk(4))])
exact_case("cdist wide", lambda mk: [(mk(1) - mk(0)) * (mk(1) - mk(0))
                                     + (mk(2) - mk(0)) * (mk(2) - mk(0))
                                     + (mk(3) - mk(0)) * (mk(3) - mk(0))])
# the case deliberately NOT used in the all-width tests: sqrt(125) is irrational
ok, _ = exact_ops_only(lambda mk: [(mk(3) * mk(3) + mk(4) * mk(4) + mk(6) * mk(6) + mk(8) * mk(8)).sqrt()])
CHECKS[0] += 1
if ok:
    FAILURES.append("fro of [[3,4],[6,8]] should NOT be exact; the test design assumes it is not")
    print("  FAIL fro [[3,4],[6,8]] unexpectedly exact")

# ======================================================== TASK: /lagoon-sort-dist

print("== lagoon-sort-dist: sorting and selection ==")


def rank_slice(xs, desc):
    """Stable order: by value, ties to the lower index."""
    return sorted(range(len(xs)), key=lambda i: (-xs[i] if desc else xs[i], i))


def sort_dim(a, dim, desc=False):
    sl, shape = slices(a, dim)
    out = [[s[i] for i in rank_slice(s, desc)] for s in sl]
    return np.moveaxis(np.array(out).reshape(shape + [-1]), -1, dim)


def argsort_dim(a, dim, desc=False):
    sl, shape = slices(a, dim)
    out = [rank_slice(s, desc) for s in sl]
    return np.moveaxis(np.array(out).reshape(shape + [-1]), -1, dim)


def argtop_dim(a, dim, k):
    sl, shape = slices(a, dim)
    out = [rank_slice(s, True)[:k] for s in sl]
    return np.moveaxis(np.array(out).reshape(shape + [k]), -1, dim)


def take_dim(a, idx, dim):
    sla, _ = slices(a, dim)
    sli, shape = slices(idx, dim)
    out = [[src[int(w)] for w in isl] for src, isl in zip(sla, sli)]
    return np.moveaxis(np.array(out).reshape(shape + [len(out[0])]), -1, dim)


S = [3, 1, 2]
T = [2, 1, 2]
G = [[3, 1], [1, 2]]
H = [[1, 5, 3], [9, 2, 7]]

check("sort S asc", flat(sort_dim(S, 0)), [1, 2, 3])
check("sort S des", flat(sort_dim(S, 0, True)), [3, 2, 1])
check("sort G dim0", flat(sort_dim(G, 0)), [1, 1, 3, 2])
check("sort G dim1", flat(sort_dim(G, 1)), [1, 3, 1, 2])
check("sort CUB dim1 des", flat(sort_dim(CUB, 1, True)), [2, 3, 0, 1, 6, 7, 4, 5])
check("argsort S asc", flat(argsort_dim(S, 0)), [1, 2, 0])
check("argsort S des", flat(argsort_dim(S, 0, True)), [0, 2, 1])
check("argsort T asc (tie -> lower index)", flat(argsort_dim(T, 0)), [1, 0, 2])
check("argsort T des (tie -> lower index)", flat(argsort_dim(T, 0, True)), [0, 2, 1])
check("argsort G dim0", flat(argsort_dim(G, 0)), [1, 0, 0, 1])
check("argtop H dim1 k2", flat(argtop_dim(H, 1, 2)), [1, 2, 0, 2])
check("argtop H dim1 k1", flat(argtop_dim(H, 1, 1)), [1, 0])
check("argtop H dim0 k1", flat(argtop_dim(H, 0, 1)), [1, 0, 1])
check("argtop H dim1 k3 == argsort des", flat(argtop_dim(H, 1, 3)), flat(argsort_dim(H, 1, True)))
check("argtop [5,5,1] k2 (tie)", flat(argtop_dim([5, 5, 1], 0, 2)), [0, 1])
check("take top-2 values", flat(take_dim(H, argtop_dim(H, 1, 2), 1)), [5, 3, 9, 7])
check("take argsort round trip", flat(take_dim(H, argsort_dim(H, 1), 1)), flat(sort_dim(H, 1)))
check("take dim0", flat(take_dim(H, [[1, 0, 0]], 0)), [9, 5, 3])
check("take repeat", flat(take_dim(H, [[1, 1, 0], [1, 1, 0]], 1)), [5, 5, 1, 2, 2, 9])
check("take CUB reverse dim1",
      flat(take_dim(CUB, [[[1, 1], [0, 0]], [[1, 1], [0, 0]]], 1)), [2, 3, 0, 1, 6, 7, 4, 5])
# numpy cross-check of the stable orders
check("argsort S asc vs numpy", flat(argsort_dim(S, 0)), flat(np.argsort(np.array(S), kind="stable")))
check("argsort T asc vs numpy", flat(argsort_dim(T, 0)), flat(np.argsort(np.array(T), kind="stable")))
check("sort G dim1 vs numpy", flat(sort_dim(G, 1)), flat(np.sort(np.array(G), axis=1)))

print("== lagoon-sort-dist: distances ==")


def cdist_sq(a, b, dt):
    a, b = np.asarray(a, dtype=dt), np.asarray(b, dtype=dt)
    out = []
    for i in range(a.shape[0]):
        row = []
        for j in range(b.shape[0]):
            t = dt(a[i, 0]) - dt(b[j, 0])
            acc = dt(t) * dt(t)
            for k in range(1, a.shape[1]):
                t = dt(a[i, k]) - dt(b[j, k])
                acc = dt(acc) + dt(dt(t) * dt(t))
            row.append(float(acc))
        out.append(row)
    return out


PA = [[0, 0], [1, 0]]
PB = [[0, 0], [0, 1], [1, 1]]
for w, dt in DT.items():
    check(f"f{w} cdist PA PB", cdist_sq(PA, PB, dt), [[0, 1, 2], [1, 2, 1]])
    check(f"f{w} cdist 3-4-5", cdist_sq([[0, 0]], [[3, 4]], dt), [[25]])
    check(f"f{w} pdist PA", cdist_sq(PA, PA, dt), [[0, 1], [1, 0]])
    check(f"f{w} pdist PB diagonal", [cdist_sq(PB, PB, dt)[i][i] for i in range(3)], [0, 0, 0])
    check(f"f{w} cdist wide", cdist_sq([[1, 2, 3]], [[0, 0, 0]], dt), [[14]])

# ==================================================== TASK: /lagoon-axis-rounding

print("== lagoon-axis-rounding: accumulation order ==")
V = [1.0, 2.0 ** -24, 2.0 ** -24]
ltr = add_f(np.float32, V)
rtl = np.float32(np.float32(V[0]) + np.float32(np.float32(V[1]) + np.float32(V[2])))
check("f32 v32 left-to-right (RNE)", bits(32, ltr), "0x3f800000")
check("f32 v32 right-to-left (RNE)", bits(32, rtl), "0x3f800001")
CHECKS[0] += 1
if bits(32, ltr) == bits(32, rtl):
    FAILURES.append("v32 fails to distinguish the fold direction")
    print("  FAIL v32 does not discriminate order")
# directed rounding, cross-checked against numpy on the nearest-even case
check("f32 v32 mode n (exact rational)", frac_bits(sum_mode(V, 32, "n"), 32), "0x3f800000")
check("f32 v32 mode z", frac_bits(sum_mode(V, 32, "z"), 32), "0x3f800000")
check("f32 v32 mode d", frac_bits(sum_mode(V, 32, "d"), 32), "0x3f800000")
check("f32 v32 mode u", frac_bits(sum_mode(V, 32, "u"), 32), "0x3f800002")
check("f32 mode n agrees with numpy", frac_bits(sum_mode(V, 32, "n"), 32), bits(32, ltr))
check("f32 ref 1.0", bits(32, 1.0), "0x3f800000")
check("f32 ref 2^-24", bits(32, 2.0 ** -24), "0x33800000")
check("f32 ref 1+2^-23", bits(32, 1 + 2.0 ** -23), "0x3f800001")
check("f32 ref 1+2^-22", bits(32, 1 + 2.0 ** -22), "0x3f800002")
# the rounder itself, against numpy, on a spread of values
probes = [0.0, 1.0, 0.5, 2.0 ** -24, 1 + 2.0 ** -23, 3.0, 1e-12, 1e12, 1 / 3, 5 / 3, 16 / 3]
for x in probes:
    check(f"rounder vs numpy f32 {x!r}", frac_bits(round_mode(Fraction(x), 32, "n"), 32), bits(32, x))
    check(f"rounder vs numpy f64 {x!r}", frac_bits(round_mode(Fraction(x), 64, "n"), 64), bits(64, x))

print("== lagoon-axis-rounding: direct form vs the Gram identity ==")


def gram_sq(a, b, dt):
    a, b = np.asarray(a, dtype=dt), np.asarray(b, dtype=dt)
    na = fold(dt, [dt(x) * dt(x) for x in a[0]], lambda p, q: p + q)
    nb = fold(dt, [dt(x) * dt(x) for x in b[0]], lambda p, q: p + q)
    dot = fold(dt, [dt(x) * dt(y) for x, y in zip(a[0], b[0])], lambda p, q: p + q)
    return float(dt(dt(na) + dt(nb)) - dt(dt(2) * dt(dot)))


GA = [[1.375, 1.875]]
GB = [[16 / 3, 5 / 3]]
for w, want_direct, want_gram in ((32, "0x417b638f", "0x417b6390"),
                                  (64, "0x402f6c71c71c71c6", "0x402f6c71c71c71c8")):
    dt = DT[w]
    d = cdist_sq(GA, GB, dt)[0][0]
    g = gram_sq(GA, GB, dt)
    check(f"f{w} cdist direct", bits(w, d), want_direct)
    check(f"f{w} cdist gram", bits(w, g), want_gram)
    CHECKS[0] += 1
    if bits(w, d) == bits(w, g):
        FAILURES.append(f"f{w} Gram case fails to discriminate the two formulas")
        print(f"  FAIL f{w} gram case does not discriminate")
    check(f"f{w} input A bits", [bits(w, x) for x in GA[0]],
          ["0x3fb00000", "0x3ff00000"] if w == 32 else
          ["0x3ff6000000000000", "0x3ffe000000000000"])
    check(f"f{w} input B bits", [bits(w, x) for x in GB[0]],
          ["0x40aaaaab", "0x3fd55555"] if w == 32 else
          ["0x4015555555555555", "0x3ffaaaaaaaaaaaab"])

print()
print(f"numpy {np.__version__}; longdouble significand bits: "
      f"{np.finfo(np.longdouble).nmant} (binary128 needs 112, so it is not usable here)")
print(f"{CHECKS[0]} checks, {len(FAILURES)} failures")
for f in FAILURES:
    print(f"  {f}")
sys.exit(1 if FAILURES else 0)
