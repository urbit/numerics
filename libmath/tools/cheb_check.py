#!/usr/bin/env python3
"""Chebyshev/minimax transcendental harness for the Hoon<->jet bit-exact effort.

Design (see numerics PR #18):
  - The algorithm uses ONLY correctly-rounded f64 primitives (+ - * /, round-to-
    int, ldexp).  For those, IEEE round-nearest-even is bit-identical in strict
    f64 (no FMA / no x87 extended precision), in Berkeley SoftFloat, and in the
    Hoon `fl` engine.  So this strict-f64 Python reference yields the exact bits
    the Hoon @rd implementation and the SoftFloat C jet must reproduce.
  - mpmath (arbitrary precision) is the INDEPENDENT truth: it designs the
    coefficients and proves the faithful-rounding (<= 1 ULP) bound.  It is NOT
    the test oracle; the algorithm-of-record (this file) is.

Run:  python3 libmath/tools/cheb_check.py exp
Requires: mpmath, numpy.
"""
import sys, struct, math
import mpmath as mp
import numpy as np
mp.mp.dps = 60

def f64(x):                       # force a Python float (strict IEEE f64)
    return struct.unpack('>d', struct.pack('>d', float(x)))[0]
def bits(x):  return struct.unpack('>Q', struct.pack('>d', f64(x)))[0]
def hexd(x):  return f"0x{bits(x):016x}"
def of_bits(b): return struct.unpack('>d', struct.pack('>Q', b))[0]

def ulps(approx, true_mpf):       # signed ULP error of f64 `approx` vs true value
    a = f64(approx)
    if a == 0: a = 0.0
    # nextafter-based ULP at a
    up = math.nextafter(a, math.inf); ulp = up - a if up != a else abs(a) * 2**-52
    if ulp == 0: ulp = 2**-1074
    return float((mp.mpf(a) - true_mpf) / mp.mpf(ulp))

# ---- exp @rd: x = k*ln2 + r, r in [-ln2/2, ln2/2]; exp = 2^k * P(r) ----
LOG2E = f64(1 / mp.log(2))
# Cody-Waite split of ln2 (fdlibm): ln2hi exact in the top bits so k*ln2hi is exact
LN2HI = f64('6.93147180369123816490e-01')
LN2LO = f64('1.90821492927058770002e-10')

def gen_exp_coeffs(deg):
    """Near-minimax monomial coeffs for exp on [-ln2/2, ln2/2] via mpmath's
    high-precision Chebyshev fit (avoids numpy lstsq ill-conditioning), each
    coefficient then rounded to f64."""
    half = mp.log(2) / 2
    cs = mp.chebyfit(lambda r: mp.e ** r, [-half, half], deg + 1)  # highest-first
    return [f64(c) for c in reversed(cs)]                          # ascending c0..cN

def horner(coeffs, r):            # ascending coeffs; strict-f64 Horner
    acc = f64(coeffs[-1])
    for c in reversed(coeffs[:-1]):
        acc = f64(f64(acc * r) + c)
    return acc

INF = math.inf
def exp_f64(x, coeffs):
    x = f64(x)
    if x != x:      return float('nan')          # NaN -> NaN
    if x == INF:    return INF                    # +inf -> +inf
    if x == -INF:   return 0.0                    # -inf -> 0
    k = int(math.floor(f64(x * LOG2E) + 0.5))     # round-nearest to int
    if k >= 1025:   return INF                     # overflow
    if k <= -1076:  return 0.0                     # underflow
    r = f64(f64(x - f64(k * LN2HI)) - f64(k * LN2LO))
    p = horner(coeffs, r)
    try:    return f64(math.ldexp(p, k))           # correctly-rounded scale (incl. subnormal)
    except OverflowError:  return INF

def check_exp():
    deg = 11
    coeffs = gen_exp_coeffs(deg)
    print(f"# exp @rd: Cody-Waite reduction + degree-{deg} minimax poly")
    print(f"LOG2E = {hexd(LOG2E)}   LN2HI = {hexd(LN2HI)}   LN2LO = {hexd(LN2LO)}")
    print("coeffs (ascending, c0..c%d), as f64 hex:" % deg)
    for i, c in enumerate(coeffs):
        print(f"  c{i:<2} = {hexd(c)}   ({c!r})")
    # faithfulness sweep
    worst = 0.0; xw = None
    xs = [mp.mpf(t)/1000 for t in range(-20000, 20001, 7)]   # x in [-20, 20]
    for xm in xs:
        x = f64(xm)
        got = exp_f64(x, coeffs)
        tru = mp.e ** mp.mpf(x)
        if tru == 0 or not math.isfinite(got): continue
        e = abs(ulps(got, tru))
        if e > worst: worst, xw = e, x
    print(f"max error over x in [-20,20]: {worst:.3f} ULP  at x={xw}")
    # expected values for the Hoon test
    print("expected (input -> output) bit patterns:")
    for x in [0.0, 0.5, 1.0, -1.0, 2.0, 10.0, -5.0, 0.1]:
        print(f"  exp({x:+}) -> {hexd(exp_f64(x, coeffs))}   in={hexd(x)}")
    print("edge cases (in -> out):")
    edges = [('+inf', INF), ('-inf', -INF), ('nan', float('nan')),
             ('709.5', 709.5), ('710.0', 710.0), ('720.0', 720.0),
             ('-744.0', -744.0), ('-745.2', -745.2), ('-750.0', -750.0)]
    for name, x in edges:
        o = exp_f64(x, coeffs)
        ib = "0x7ff0000000000000" if x==INF else "0xfff0000000000000" if x==-INF else \
             "0x7ff8000000000000" if x!=x else hexd(x)
        print(f"  exp({name:>7}) -> {hexd(o):>18}   in={ib}")

# ============================ @rs (f32) ============================
def f32(x):   return struct.unpack('>f', struct.pack('>f', float(x)))[0]
def bits32(x):return struct.unpack('>I', struct.pack('>f', f32(x)))[0]
def hexs(x):  return f"0x{bits32(x):08x}"

def ulps32(approx, true_mpf):
    a = f32(approx)
    up = f32(math.nextafter(a, math.inf)); ulp = up - a if up != a else abs(a) * 2**-23
    if ulp == 0: ulp = 2**-149
    return float((mp.mpf(a) - true_mpf) / mp.mpf(ulp))

# Cody-Waite split of ln2 for f32 (fdlibm expf): low mantissa bits of HI are 0
LOG2E_S = f32(1 / mp.log(2))
LN2HI_S = struct.unpack('>f', struct.pack('>I', 0x3f317200))[0]  # 6.9314575195e-1
LN2LO_S = struct.unpack('>f', struct.pack('>I', 0x35bfbe8e))[0]  # 1.4286067653e-6

def gen_exp_coeffs_f32(deg):
    half = mp.log(2) / 2
    cs = mp.chebyfit(lambda r: mp.e ** r, [-half, half], deg + 1)
    return [f32(c) for c in reversed(cs)]

def horner32(coeffs, r):
    acc = f32(coeffs[-1])
    for c in reversed(coeffs[:-1]):
        acc = f32(f32(acc * r) + c)
    return acc

def exp_f32(x, coeffs):
    x = f32(x)
    if x != x:      return float('nan')
    if x == INF:    return INF
    if x == -INF:   return 0.0
    k = int(math.floor(f32(x * LOG2E_S) + 0.5))
    if k >= 129:    return INF                     # overflow (f32 max exp ~88.7)
    if k <= -151:   return 0.0                     # underflow (smallest subnormal ~2^-149)
    r = f32(f32(x - f32(k * LN2HI_S)) - f32(k * LN2LO_S))
    p = horner32(coeffs, r)
    try:    return f32(math.ldexp(p, k))
    except OverflowError:  return INF

def check_exp_rs():
    deg = 6
    coeffs = gen_exp_coeffs_f32(deg)
    print(f"# exp @rs: Cody-Waite reduction + degree-{deg} minimax poly (f32)")
    print(f"LOG2E = {hexs(LOG2E_S)}   LN2HI = {hexs(LN2HI_S)}   LN2LO = {hexs(LN2LO_S)}")
    print("coeffs (ascending, c0..c%d), as f32 hex:" % deg)
    for i, c in enumerate(coeffs):
        print(f"  c{i:<2} = {hexs(c)}   ({c!r})")
    worst = 0.0; xw = None
    for t in range(-20000, 20001, 7):
        x = f32(mp.mpf(t) / 1000)
        got = exp_f32(x, coeffs); tru = mp.e ** mp.mpf(x)
        if tru == 0 or not math.isfinite(got): continue
        e = abs(ulps32(got, tru))
        if e > worst: worst, xw = e, x
    print(f"max error over x in [-20,20]: {worst:.3f} ULP  at x={xw}")
    print("expected (input -> output) bit patterns:")
    for x in [0.0, 0.5, 1.0, -1.0, 2.0, 10.0, -5.0, 0.1]:
        print(f"  exp({x:+}) -> {hexs(exp_f32(x, coeffs))}   in={hexs(x)}")
    print("edge cases (in -> out):")
    edges = [('+inf', INF), ('-inf', -INF), ('nan', float('nan')),
             ('88.0', 88.0), ('89.0', 89.0), ('100.0', 100.0),
             ('-103.0', -103.0), ('-104.0', -104.0), ('-110.0', -110.0)]
    for name, x in edges:
        o = exp_f32(x, coeffs)
        ib = "0x7f800000" if x==INF else "0xff800000" if x==-INF else \
             "0x7fc00000" if x!=x else hexs(x)
        print(f"  exp({name:>7}) -> {hexs(o):>12}   in={ib}")

# ============================ log @rd ============================
#  x = 2^e * m, m in [sqrt(1/2), sqrt(2)); log = e*ln2 + 2s*P(z),
#  s = (m-1)/(m+1), z = s*s, P(z) = sum_k z^k/(2k+1) (atanh series, z<=0.0294).
SQRT2  = f64(mp.sqrt(2))
ONE    = f64(1.0)
def gen_log_coeffs(deg):                       # P2(z) = 1/3 + z/5 + z^2/7 + ...
    return [f64(mp.mpf(1) / (2*k + 3)) for k in range(deg + 1)]

def log_f64(x, coeffs):
    x = f64(x)
    if x != x:                 return float('nan')        # NaN
    if x == INF:               return INF                  # +inf
    if x < 0.0 or x == -INF:   return float('nan')         # x<0 -> NaN
    if x == 0.0:               return -INF                 # log(+-0) -> -inf
    # normalise subnormals so the bit-extraction reduction sees a normal m
    add_e = 0
    if x < 2.2250738585072014e-308:                        # < smallest normal
        x = f64(x * 18014398509481984.0); add_e = -54      # *2^54
    b = bits(x); ef = ((b >> 52) & 0x7ff) - 1023; m = of_bits((b & ((1<<52)-1)) | 0x3ff0000000000000)
    if m >= SQRT2:  m = f64(m * 0.5); ef += 1
    ef += add_e
    # log(1+f) = f - s*(f - R), R = 2z*P2(z): keeps f as the exact leading term
    f = f64(m - ONE)
    s = f64(f / f64(m + ONE))
    z = f64(s * s)
    p2 = horner(coeffs, z)
    r  = f64(f64(f64(z + z)) * p2)                         # 2z*P2(z)
    l1 = f64(f - f64(s * f64(f - r)))                      # log(1+f)
    e = f64(float(ef))
    return f64(f64(e * LN2HI) + f64(l1 + f64(e * LN2LO)))

def check_log():
    deg = 9
    coeffs = gen_log_coeffs(deg)
    print(f"# log @rd: x=2^e*m reduction + degree-{deg} atanh poly")
    print(f"LN2HI = {hexd(LN2HI)}   LN2LO = {hexd(LN2LO)}   SQRT2 = {hexd(SQRT2)}")
    print("coeffs (ascending, c0..c%d), as f64 hex:" % deg)
    for i, c in enumerate(coeffs):
        print(f"  c{i:<2} = {hexd(c)}   ({c!r})")
    worst = 0.0; xw = None
    for t in range(1, 200001, 7):
        x = f64(mp.mpf(t) / 1000)
        got = log_f64(x, coeffs); tru = mp.log(mp.mpf(x))
        if tru == 0 or not math.isfinite(got): continue
        e = abs(ulps(got, tru))
        if e > worst: worst, xw = e, x
    print(f"max error over x in (0,200]: {worst:.3f} ULP  at x={xw}")
    print("expected (input -> output) bit patterns:")
    for x in [1.0, 2.0, 0.5, 10.0, 100.0, 0.1, 1.0e-300, 7.389056098930650]:
        print(f"  log({x}) -> {hexd(log_f64(x, coeffs))}   in={hexd(x)}")
    print("edge cases (in -> out):")
    edges = [('+inf', INF), ('-inf', -INF), ('nan', float('nan')),
             ('0.0', 0.0), ('-1.0', -1.0), ('1.0', 1.0)]
    for name, x in edges:
        o = log_f64(x, coeffs)
        ib = "0x7ff0000000000000" if x==INF else "0xfff0000000000000" if x==-INF else \
             "0x7ff8000000000000" if x!=x else hexd(x)
        print(f"  log({name:>6}) -> {hexd(o):>18}   in={ib}")

# ============================ log @rs ============================
SQRT2_S = f32(mp.sqrt(2))
def gen_log_coeffs_f32(deg):
    return [f32(mp.mpf(1) / (2*k + 3)) for k in range(deg + 1)]

def log_f32(x, coeffs):
    x = f32(x)
    if x != x:                 return float('nan')
    if x == INF:               return INF
    if x < 0.0 or x == -INF:   return float('nan')
    if x == 0.0:               return -INF
    add_e = 0
    if x < 1.1754943508222875e-38:                         # < smallest normal f32
        x = f32(x * 16777216.0); add_e = -24               # *2^24
    b = bits32(x); ef = ((b >> 23) & 0xff) - 127
    m = struct.unpack('>f', struct.pack('>I', (b & 0x7fffff) | 0x3f800000))[0]
    if m >= SQRT2_S:  m = f32(m * 0.5); ef += 1
    ef += add_e
    f = f32(m - 1.0)
    s = f32(f / f32(m + 1.0))
    z = f32(s * s)
    p2 = horner32(coeffs, z)
    r  = f32(f32(z + z) * p2)
    l1 = f32(f - f32(s * f32(f - r)))
    e = f32(float(ef))
    return f32(f32(e * LN2HI_S) + f32(l1 + f32(e * LN2LO_S)))

def check_log_rs():
    deg = 4
    coeffs = gen_log_coeffs_f32(deg)
    print(f"# log @rs: x=2^e*m reduction + degree-{deg} atanh poly (f32)")
    print(f"LN2HI = {hexs(LN2HI_S)}   LN2LO = {hexs(LN2LO_S)}   SQRT2 = {hexs(SQRT2_S)}")
    print("coeffs (ascending, c0..c%d), as f32 hex:" % deg)
    for i, c in enumerate(coeffs):
        print(f"  c{i:<2} = {hexs(c)}   ({c!r})")
    worst = 0.0; xw = None
    for t in range(1, 200001, 7):
        x = f32(mp.mpf(t) / 1000)
        got = log_f32(x, coeffs); tru = mp.log(mp.mpf(x))
        if tru == 0 or not math.isfinite(got): continue
        e = abs(ulps32(got, tru))
        if e > worst: worst, xw = e, x
    print(f"max error over x in (0,200]: {worst:.3f} ULP  at x={xw}")
    print("expected (input -> output) bit patterns:")
    for x in [1.0, 2.0, 0.5, 10.0, 100.0, 0.1, 1.0e-40, 7.389056]:
        print(f"  log({x}) -> {hexs(log_f32(x, coeffs))}   in={hexs(x)}")
    print("edge cases (in -> out):")
    edges = [('+inf', INF), ('-inf', -INF), ('nan', float('nan')),
             ('0.0', 0.0), ('-1.0', -1.0), ('1.0', 1.0)]
    for name, x in edges:
        o = log_f32(x, coeffs)
        ib = "0x7f800000" if x==INF else "0xff800000" if x==-INF else \
             "0x7fc00000" if x!=x else hexs(x)
        print(f"  log({name:>6}) -> {hexs(o):>12}   in={ib}")

# ============================ sin/cos @rd ============================
#  x = q*(pi/2) + r, r in [-pi/4, pi/4]; pick sin/cos kernel by q&3.
PIO2    = mp.pi / 2
INVPIO2 = f64(2 / mp.pi)
PIO2_1  = of_bits(bits(f64(PIO2)) & ~((1 << 22) - 1))   # pi/2, low 22 mantissa bits 0
PIO2_1T = f64(PIO2 - mp.mpf(PIO2_1))                    # tail
SIN_C = [f64(mp.mpf((-1)**(k+1)) / mp.factorial(2*k+3)) for k in range(8)]  # S(z): -1/6,1/120,..
COS_C = [f64(mp.mpf((-1)**k)     / mp.factorial(2*k+4)) for k in range(8)]  # C(z):  1/24,-1/720,..

def ksin(x, y):                               # fdlibm __kernel_sin(x,y); x+y = reduced arg
    z = f64(x*x)
    r = horner(SIN_C[1:], z)                   # S2 + z*S3 + ...  (exact-Taylor, 7 terms)
    v = f64(z*x)
    return f64(x - f64(f64(f64(z*f64(f64(0.5*y) - f64(v*r))) - y) - f64(v*SIN_C[0])))
def kcos(x, y):                               # fdlibm __kernel_cos(x,y)
    z = f64(x*x)
    rc = horner(COS_C, z)                      # C1 + z*C2 + ...  (z^2*rc = r^4/24 - ...)
    hz = f64(0.5*z); w2 = f64(1.0 - hz)
    return f64(w2 + f64(f64(f64(1.0 - w2) - hz) + f64(f64(f64(z*z)*rc) - f64(x*y))))
def reduce_pio2(x):                           # x = q*pi/2 + (rhi+rlo); Fast2Sum tail
    q = int(math.floor(f64(x * INVPIO2) + 0.5))
    t = f64(x - f64(q * PIO2_1))              # exact in the Sterbenz region
    w = f64(q * PIO2_1T)
    rhi = f64(t - w)
    rlo = f64(f64(t - rhi) - w)
    return q, rhi, rlo
def sin_f64(x):                                # compute on |x|, apply odd symmetry
    x = f64(x)
    if x != x or x == INF or x == -INF: return float('nan')
    if x == 0.0: return x                      # +-0 -> +-0
    neg = bits(x) >> 63; ax = f64(abs(x))
    q, rhi, rlo = reduce_pio2(ax); m = q & 3
    v = [ksin(rhi, rlo), kcos(rhi, rlo),
         f64(-ksin(rhi, rlo)), f64(-kcos(rhi, rlo))][m]
    return f64(-v) if neg else v
def cos_f64(x):                                # cos is even
    x = f64(x)
    if x != x or x == INF or x == -INF: return float('nan')
    ax = f64(abs(x))
    q, rhi, rlo = reduce_pio2(ax); m = q & 3
    return [kcos(rhi, rlo), f64(-ksin(rhi, rlo)),
            f64(-kcos(rhi, rlo)), ksin(rhi, rlo)][m]

# ---- tan @rd: fdlibm __kernel_tan over the q*pi/2 reduction ----
TAN_T = [f64(s) for s in [
    '3.33333333333334091986e-01','1.33333333333201242699e-01',
    '5.39682539762260521377e-02','2.18694882948595424599e-02',
    '8.86323982359930005737e-03','3.59207910759131235356e-03',
    '1.45620945432529025516e-03','5.88041240820264096874e-04',
    '2.46463134818469906812e-04','7.81794442939557092300e-05',
    '7.14072491382608190305e-05','-1.85586374855275456654e-05',
    '2.59073051863633712884e-05']]
PIO4 = f64(mp.pi/4); PIO4LO = of_bits(0x3C81A62633145C07)
TAN_BIG = of_bits(0x3FE5942800000000)          # ~0.6744 (fdlibm high-word cut)
def head0(x): return of_bits(bits(x) & 0xffffffff00000000)
def ktan(x, y, iy):
    hx_neg = bits(x) >> 63
    big = f64(abs(x)) >= TAN_BIG
    if big:
        if hx_neg: x = f64(-x); y = f64(-y)
        z = f64(PIO4 - x); w0 = f64(PIO4LO - y); x = f64(z + w0); y = 0.0
    z = f64(x*x); w = f64(z*z)
    r = horner([TAN_T[1],TAN_T[3],TAN_T[5],TAN_T[7],TAN_T[9],TAN_T[11]], w)
    v = f64(z * horner([TAN_T[2],TAN_T[4],TAN_T[6],TAN_T[8],TAN_T[10],TAN_T[12]], w))
    s = f64(z * x)
    r = f64(y + f64(z * f64(f64(s * f64(r + v)) + y)))
    r = f64(r + f64(TAN_T[0] * s))
    w = f64(x + r)
    if big:
        fac = 1.0 if not hx_neg else -1.0
        v = float(iy)
        return f64(fac * f64(v - f64(2.0 * f64(x - f64(f64(f64(w*w) / f64(w+v)) - r)))))
    if iy == 1: return w
    z2 = head0(w); v2 = f64(r - f64(z2 - x)); a = f64(-1.0 / w); t2 = head0(a)
    s2 = f64(1.0 + f64(t2 * z2))
    return f64(t2 + f64(a * f64(s2 + f64(t2 * v2))))
def tan_f64(x):
    x = f64(x)
    if x != x or x == INF or x == -INF: return float('nan')
    if x == 0.0: return x
    neg = bits(x) >> 63; ax = f64(abs(x))
    q, rhi, rlo = reduce_pio2(ax)
    iy = -1 if (q & 1) else 1
    t = ktan(rhi, rlo, iy)
    return f64(-t) if neg else t

def check_tan():
    print("# tan @rd: q*pi/2 reduction + fdlibm __kernel_tan (deg-13 + cot path)")
    print("T: " + " ".join(hexd(c) for c in TAN_T))
    print(f"PIO4={hexd(PIO4)} PIO4LO={hexd(PIO4LO)} BIG={hexd(TAN_BIG)}")
    wt = wr = 0.0; xw = None
    for t in range(-200000, 200001, 7):
        x = f64(mp.mpf(t) / 1000)
        tr = mp.tan(mp.mpf(x))
        if not math.isfinite(float(tr)) or abs(tr) > 1e15 or abs(tr) < 1e-9: continue
        gt = tan_f64(x); gr = f64(sin_f64(x) / cos_f64(x))   # dedicated vs ratio
        if math.isfinite(gt):
            e = abs(ulps(gt, tr))
            if e > wt: wt, xw = e, x
        if math.isfinite(gr): wr = max(wr, abs(ulps(gr, tr)))
    print(f"dedicated max {wt:.3f} ULP at x={xw};  sin/cos ratio max {wr:.3f} ULP")
    for x in [0.0, 0.5, 1.0, -1.0, 0.7853981633974483, 2.0, 10.0, 100.0]:
        print(f"  tan({x}) -> {hexd(tan_f64(x))}   in={hexd(x)}")
    for name, x in [('+inf', INF), ('nan', float('nan')), ('-0', -0.0)]:
        o = tan_f64(x)
        ib = "0x7ff0000000000000" if x==INF else "0x7ff8000000000000" if x!=x else hexd(x)
        print(f"  tan({name}) -> {hexd(o)}   in={ib}")

def check_trig(which):
    fn = sin_f64 if which == 'sin' else cos_f64
    tru = mp.sin if which == 'sin' else mp.cos
    print(f"# {which} @rd: x=q*pi/2+r reduction; kernels deg-7 (sin)/deg-7 (cos)")
    print(f"INVPIO2={hexd(INVPIO2)}  PIO2_1={hexd(PIO2_1)}  PIO2_1T={hexd(PIO2_1T)}")
    nm = 'SIN_C' if which=='sin' else 'COS_C'  # both kernels are always needed; print both
    for label, arr in [('SIN_C', SIN_C), ('COS_C', COS_C)]:
        print(f"{label}: " + " ".join(hexd(c) for c in arr))
    worst = 0.0; xw = None
    for t in range(-200000, 200001, 7):
        x = f64(mp.mpf(t) / 1000)
        got = fn(x); tr = tru(mp.mpf(x))
        if not math.isfinite(got) or abs(tr) < 1e-9: continue
        e = abs(ulps(got, tr))
        if e > worst: worst, xw = e, x
    print(f"max error over x in [-200,200] (|f|>1e-9): {worst:.3f} ULP  at x={xw}")
    print("expected (input -> output):")
    for x in [0.0, 0.5, 1.0, -1.0, 1.5707963267948966, 3.141592653589793, 10.0, 100.0]:
        print(f"  {which}({x}) -> {hexd(fn(x))}   in={hexd(x)}")
    print("edges:")
    for name, x in [('+inf', INF), ('-inf', -INF), ('nan', float('nan')), ('-0', -0.0)]:
        o = fn(x)
        ib = "0x7ff0000000000000" if x==INF else "0xfff0000000000000" if x==-INF else \
             "0x7ff8000000000000" if x!=x else hexd(x)
        print(f"  {which}({name:>5}) -> {hexd(o):>18}   in={ib}")

# ============================ sin/cos @rs ============================
def of_bits32(b): return struct.unpack('>f', struct.pack('>I', b))[0]
INVPIO2_S = f32(2 / mp.pi)
# 3-part pi/2: P1,P2 each ~10 sig bits (so q*Pi exact for q<2^14), P3 the tail.
PIO2_1_S  = of_bits32(bits32(f32(PIO2)) & ~0x1fff)
PIO2_2_S  = of_bits32(bits32(f32(mp.mpf(PIO2) - mp.mpf(PIO2_1_S))) & ~0x1fff)
PIO2_3_S  = f32(mp.mpf(PIO2) - mp.mpf(PIO2_1_S) - mp.mpf(PIO2_2_S))
SIN_C_S = [f32(mp.mpf((-1)**(k+1)) / mp.factorial(2*k+3)) for k in range(5)]
COS_C_S = [f32(mp.mpf((-1)**k)     / mp.factorial(2*k+4)) for k in range(5)]

def ksin32(x, y):
    z = f32(x*x); r = horner32(SIN_C_S[1:], z); v = f32(z*x)
    return f32(x - f32(f32(f32(z*f32(f32(0.5*y) - f32(v*r))) - y) - f32(v*SIN_C_S[0])))
def kcos32(x, y):
    z = f32(x*x); rc = horner32(COS_C_S, z)
    hz = f32(0.5*z); w2 = f32(1.0 - hz)
    return f32(w2 + f32(f32(f32(1.0 - w2) - hz) + f32(f32(f32(z*z)*rc) - f32(x*y))))
def reduce_pio2_32(x):
    q = int(math.floor(f32(x * INVPIO2_S) + 0.5))
    r = f32(x - f32(q * PIO2_1_S)); r = f32(r - f32(q * PIO2_2_S))
    w = f32(q * PIO2_3_S); rhi = f32(r - w); rlo = f32(f32(r - rhi) - w)
    return q, rhi, rlo
def sin_f32(x):
    x = f32(x)
    if x != x or x == INF or x == -INF: return float('nan')
    if x == 0.0: return x
    neg = bits32(x) >> 31; ax = f32(abs(x))
    q, rhi, rlo = reduce_pio2_32(ax); m = q & 3
    v = [ksin32(rhi,rlo), kcos32(rhi,rlo), f32(-ksin32(rhi,rlo)), f32(-kcos32(rhi,rlo))][m]
    return f32(-v) if neg else v
def cos_f32(x):
    x = f32(x)
    if x != x or x == INF or x == -INF: return float('nan')
    ax = f32(abs(x))
    q, rhi, rlo = reduce_pio2_32(ax); m = q & 3
    return [kcos32(rhi,rlo), f32(-ksin32(rhi,rlo)), f32(-kcos32(rhi,rlo)), ksin32(rhi,rlo)][m]

def check_trig_rs(which):
    fn = sin_f32 if which == 'sin' else cos_f32
    tru = mp.sin if which == 'sin' else mp.cos
    print(f"# {which} @rs: x=q*pi/2+r reduction (f32)")
    print(f"INVPIO2={hexs(INVPIO2_S)}  PIO2_1={hexs(PIO2_1_S)}  PIO2_2={hexs(PIO2_2_S)}  PIO2_3={hexs(PIO2_3_S)}")
    print("SIN_C: " + " ".join(hexs(c) for c in SIN_C_S))
    print("COS_C: " + " ".join(hexs(c) for c in COS_C_S))
    worst = 0.0; xw = None
    for t in range(-200000, 200001, 7):
        x = f32(mp.mpf(t) / 1000); got = fn(x); tr = tru(mp.mpf(x))
        if not math.isfinite(got) or abs(tr) < 1e-6: continue
        e = abs(ulps32(got, tr))
        if e > worst: worst, xw = e, x
    print(f"max error over x in [-200,200] (|f|>1e-6): {worst:.3f} ULP  at x={xw}")
    print("expected (input -> output):")
    for x in [0.0, 0.5, 1.0, -1.0, 1.5707964, 3.1415927, 10.0, 100.0]:
        print(f"  {which}({x}) -> {hexs(fn(x))}   in={hexs(x)}")
    for name, x in [('+inf', INF), ('-inf', -INF), ('nan', float('nan')), ('-0', -0.0)]:
        o = fn(x)
        ib = "0x7f800000" if x==INF else "0xff800000" if x==-INF else \
             "0x7fc00000" if x!=x else hexs(x)
        print(f"  {which}({name:>5}) -> {hexs(o):>12}   in={ib}")

# ============================ atan @rd ============================
#  fdlibm s_atan: reduce |x| against breakpoints 7/16,11/16,19/16,39/16 to a
#  small argument near atan(0.5)/atan(1)/atan(1.5)/atan(inf), minimax poly.
ATAN_AT = [f64(s) for s in [
    '3.33333333333329318027e-01','-1.99999999998764832476e-01',
    '1.42857142725034663711e-01','-1.11111104054623557880e-01',
    '9.09088713343650656196e-02','-7.69187620504482999495e-02',
    '6.66107313738753120669e-02','-5.83357013379057348645e-02',
    '4.97687799461593236017e-02','-3.65315727442169155270e-02',
    '1.62858201153657823623e-02']]
ATAN_BP = [f64(7)/16, f64(11)/16, f64(19)/16, f64(39)/16]
ATANHI = [f64(mp.atan(mp.mpf('0.5'))), f64(mp.pi/4), f64(mp.atan(mp.mpf('1.5'))), f64(mp.pi/2)]
ATANLO = [f64(mp.atan(mp.mpf('0.5')) - mp.mpf(ATANHI[0])), f64(mp.pi/4 - mp.mpf(ATANHI[1])),
          f64(mp.atan(mp.mpf('1.5')) - mp.mpf(ATANHI[2])), f64(mp.pi/2 - mp.mpf(ATANHI[3]))]

def atan_kernel(x):                            # x >= 0 finite
    if x < ATAN_BP[0]:      idd = -1; xr = x
    elif x < ATAN_BP[1]:    idd = 0;  xr = f64(f64(f64(x+x)-ONE) / f64(2.0+x))
    elif x < ATAN_BP[2]:    idd = 1;  xr = f64(f64(x-ONE) / f64(x+ONE))
    elif x < ATAN_BP[3]:    idd = 2;  xr = f64(f64(x-1.5) / f64(ONE+f64(1.5*x)))
    else:                   idd = 3;  xr = f64(-1.0 / x)
    z = f64(xr*xr)
    s = f64(z * horner(ATAN_AT, z))            # (s1+s2)
    if idd < 0:  return f64(xr - f64(xr*s))
    return f64(ATANHI[idd] - f64(f64(f64(xr*s) - ATANLO[idd]) - xr))
def atan_f64(x):
    x = f64(x)
    if x != x:    return float('nan')
    if x == INF:  return ATANHI[3]
    if x == -INF: return f64(-ATANHI[3])
    if x == 0.0:  return x
    neg = bits(x) >> 63; ax = f64(abs(x))
    r = atan_kernel(ax)
    return f64(-r) if neg else r

def check_atan():
    print("# atan @rd: fdlibm breakpoint reduction + degree-10 minimax poly")
    print("AT: " + " ".join(hexd(c) for c in ATAN_AT))
    print("BP: " + " ".join(hexd(c) for c in ATAN_BP))
    print("ATANHI: " + " ".join(hexd(c) for c in ATANHI))
    print("ATANLO: " + " ".join(hexd(c) for c in ATANLO))
    worst = 0.0; xw = None
    for t in range(-500000, 500001, 7):
        x = f64(mp.mpf(t) / 1000); got = atan_f64(x); tr = mp.atan(mp.mpf(x))
        if abs(tr) < 1e-12 or not math.isfinite(got): continue
        e = abs(ulps(got, tr))
        if e > worst: worst, xw = e, x
    print(f"max error over x in [-500,500]: {worst:.3f} ULP  at x={xw}")
    print("expected:")
    for x in [0.0, 0.5, 1.0, -1.0, 1.5, 2.0, 10.0, 0.1, -0.7]:
        print(f"  atan({x}) -> {hexd(atan_f64(x))}   in={hexd(x)}")
    for name, x in [('+inf', INF), ('-inf', -INF), ('nan', float('nan')), ('-0', -0.0)]:
        o = atan_f64(x)
        ib = "0x7ff0000000000000" if x==INF else "0xfff0000000000000" if x==-INF else \
             "0x7ff8000000000000" if x!=x else hexd(x)
        print(f"  atan({name:>5}) -> {hexd(o):>18}   in={ib}")

# ============================ atan @rs ============================
ATAN_AT_S = [f32(s) for s in ['3.3333328366e-01','-1.9999158382e-01',
    '1.4253635705e-01','-1.0648017377e-01','6.1687607318e-02']]
ATAN_BP_S = [f32(7)/16, f32(11)/16, f32(19)/16, f32(39)/16]
ATANHI_S = [f32(mp.atan(mp.mpf('0.5'))), f32(mp.pi/4), f32(mp.atan(mp.mpf('1.5'))), f32(mp.pi/2)]
ATANLO_S = [f32(mp.atan(mp.mpf('0.5')) - mp.mpf(ATANHI_S[0])), f32(mp.pi/4 - mp.mpf(ATANHI_S[1])),
            f32(mp.atan(mp.mpf('1.5')) - mp.mpf(ATANHI_S[2])), f32(mp.pi/2 - mp.mpf(ATANHI_S[3]))]
def atan_kernel32(x):
    if x < ATAN_BP_S[0]:   idd = -1; xr = x
    elif x < ATAN_BP_S[1]: idd = 0;  xr = f32(f32(f32(x+x)-1.0) / f32(2.0+x))
    elif x < ATAN_BP_S[2]: idd = 1;  xr = f32(f32(x-1.0) / f32(x+1.0))
    elif x < ATAN_BP_S[3]: idd = 2;  xr = f32(f32(x-1.5) / f32(1.0+f32(1.5*x)))
    else:                  idd = 3;  xr = f32(-1.0 / x)
    z = f32(xr*xr)
    s = f32(z * horner32(ATAN_AT_S, z))
    if idd < 0:  return f32(xr - f32(xr*s))
    return f32(ATANHI_S[idd] - f32(f32(f32(xr*s) - ATANLO_S[idd]) - xr))
def atan_f32(x):
    x = f32(x)
    if x != x:    return float('nan')
    if x == INF:  return ATANHI_S[3]
    if x == -INF: return f32(-ATANHI_S[3])
    if x == 0.0:  return x
    neg = bits32(x) >> 31; ax = f32(abs(x))
    r = atan_kernel32(ax)
    return f32(-r) if neg else r

def check_atan_rs():
    print("# atan @rs: fdlibm breakpoint reduction + degree-4 minimax poly (f32)")
    print("AT: " + " ".join(hexs(c) for c in ATAN_AT_S))
    print("BP: " + " ".join(hexs(c) for c in ATAN_BP_S))
    print("ATANHI: " + " ".join(hexs(c) for c in ATANHI_S))
    print("ATANLO: " + " ".join(hexs(c) for c in ATANLO_S))
    worst = 0.0; xw = None
    for t in range(-500000, 500001, 7):
        x = f32(mp.mpf(t) / 1000); got = atan_f32(x); tr = mp.atan(mp.mpf(x))
        if abs(tr) < 1e-7 or not math.isfinite(got): continue
        e = abs(ulps32(got, tr))
        if e > worst: worst, xw = e, x
    print(f"max error over x in [-500,500]: {worst:.3f} ULP  at x={xw}")
    print("expected:")
    for x in [0.0, 0.5, 1.0, -1.0, 1.5, 2.0, 10.0, 0.1, -0.7]:
        print(f"  atan({x}) -> {hexs(atan_f32(x))}   in={hexs(x)}")
    for name, x in [('+inf', INF), ('-inf', -INF), ('nan', float('nan')), ('-0', -0.0)]:
        o = atan_f32(x)
        ib = "0x7f800000" if x==INF else "0xff800000" if x==-INF else \
             "0x7fc00000" if x!=x else hexs(x)
        print(f"  atan({name:>5}) -> {hexs(o):>12}   in={ib}")

# ============================ asin @rd ============================
#  fdlibm e_asin: |x|<0.5 rational x+x*R(x^2); |x| in [0.5,1) via t=(1-|x|)/2,
#  s=sqrt(t), asin = pi/2 - 2*(s + s*R(t)).  R = P/Q (pS/qS coeffs).
PS = [f64(s) for s in ['1.66666666666666657415e-01','-3.25565818622400915405e-01',
    '2.01212532134862925881e-01','-4.00555345006794114027e-02',
    '7.91534994289814532176e-04','3.47933107596021167570e-05']]
QS = [f64(s) for s in ['-2.40339491173441421878e+00','2.02094576023350569471e+00',
    '-6.88283971605453293030e-01','7.70381505559019352791e-02']]
PIO2_H = f64(mp.pi/2); PIO2_L = f64(mp.pi/2 - mp.mpf(PIO2_H)); PIO4_H = f64(mp.pi/4)
ASIN_THRESH = of_bits(0x3fef333300000000)      # ~0.975 (fdlibm high-word cut)

def asin_R(t):                                 # P(t)/Q(t)
    p = f64(t * horner(PS, t))
    q = f64(ONE + f64(t * horner(QS, t)))
    return f64(p / q)
def asin_f64(x):
    x = f64(x)
    if x != x: return float('nan')
    ax = f64(abs(x)); sgn = bits(x) >> 63
    if ax > ONE: return float('nan')
    if ax == ONE:
        return f64(f64(x * PIO2_H) + f64(x * PIO2_L))    # +-pi/2 with sign of x
    if ax < 0.5:
        if ax < 1.49e-8: return x                        # |x|<2^-26: asin(x)=x
        t = f64(x * x)
        return f64(x + f64(x * asin_R(t)))
    w = f64(ONE - ax); t = f64(w * 0.5)
    r = asin_R(t); s = f64(math.sqrt(t))
    if ax >= ASIN_THRESH:                                # near 1: simple form
        res = f64(PIO2_H - f64(f64(2.0 * f64(s + f64(s * r))) - PIO2_L))
    else:                                                # head/tail recovers low bits of s
        df = of_bits(bits(s) & 0xffffffff00000000)
        c  = f64(f64(t - f64(df * df)) / f64(s + df))
        p2 = f64(f64(2.0 * f64(s * r)) - f64(PIO2_L - f64(2.0 * c)))
        q2 = f64(PIO4_H - f64(2.0 * df))
        res = f64(PIO4_H - f64(p2 - q2))
    return f64(-res) if sgn else res

PI_H = f64(mp.pi)
def acos_f64(x):                               # fdlibm e_acos
    x = f64(x)
    if x != x: return float('nan')
    ax = f64(abs(x)); neg = bits(x) >> 63
    if ax > ONE: return float('nan')
    if ax == ONE:
        return 0.0 if not neg else f64(PI_H + f64(2.0 * PIO2_L))   # acos(1)=0, acos(-1)=pi
    if ax < 0.5:
        if ax < 6.94e-18: return PIO2_H                            # |x|<2^-57
        z = f64(x * x); r = asin_R(z)
        return f64(PIO2_H - f64(x - f64(PIO2_L - f64(x * r))))
    if neg:                                                        # x <= -0.5
        z = f64(f64(ONE + x) * 0.5); s = f64(math.sqrt(z)); r = asin_R(z)
        w = f64(f64(r * s) - PIO2_L)
        return f64(PI_H - f64(2.0 * f64(s + w)))
    z = f64(f64(ONE - x) * 0.5); s = f64(math.sqrt(z))             # x >= 0.5
    df = of_bits(bits(s) & 0xffffffff00000000)
    c  = f64(f64(z - f64(df * df)) / f64(s + df))
    r  = asin_R(z); w = f64(f64(r * s) + c)
    return f64(2.0 * f64(df + w))

def check_acos():
    print("# acos @rd: fdlibm rational kernel (shares asin P/Q)")
    worst = 0.0; xw = None
    for t in range(-1000000, 1000001, 3):
        x = f64(mp.mpf(t) / 1000000)
        got = acos_f64(x); tr = mp.acos(mp.mpf(x))
        if abs(tr) < 1e-12 or not math.isfinite(got): continue
        e = abs(ulps(got, tr))
        if e > worst: worst, xw = e, x
    print(f"max error over x in [-1,1]: {worst:.3f} ULP  at x={xw}")
    for x in [0.0, 0.5, 1.0, -1.0, 0.25, 0.75, 0.9, -0.6, -0.9, 0.1]:
        print(f"  acos({x}) -> {hexd(acos_f64(x))}   in={hexd(x)}")
    for name, x in [('nan', float('nan')), ('1.5', 1.5), ('-2', -2.0)]:
        print(f"  acos({name}) -> {hexd(acos_f64(x))}")

def check_asin():
    print("# asin @rd: fdlibm rational kernel")
    print("PS: " + " ".join(hexd(c) for c in PS))
    print("QS: " + " ".join(hexd(c) for c in QS))
    print(f"PIO2_H={hexd(PIO2_H)} PIO2_L={hexd(PIO2_L)}")
    worst = 0.0; xw = None
    for t in range(-1000000, 1000001, 3):
        x = f64(mp.mpf(t) / 1000000)
        got = asin_f64(x); tr = mp.asin(mp.mpf(x))
        if abs(tr) < 1e-12 or not math.isfinite(got): continue
        e = abs(ulps(got, tr))
        if e > worst: worst, xw = e, x
    print(f"max error over x in [-1,1]: {worst:.3f} ULP  at x={xw}")
    for x in [0.0, 0.5, 1.0, -1.0, 0.25, 0.75, 0.9, 0.99, -0.6, 0.1]:
        print(f"  asin({x}) -> {hexd(asin_f64(x))}   in={hexd(x)}")
    for name, x in [('nan', float('nan')), ('-0', -0.0), ('1.5', 1.5), ('-2', -2.0)]:
        print(f"  asin({name}) -> {hexd(asin_f64(x))}")

# ============================ asin/acos @rs ============================
PS_S = [f32(s) for s in ['1.6666586697e-01','-4.2743422091e-02','-8.6563630030e-03']]
QS1_S = f32('-7.0662963390e-01')
PIO2_HS = f32(mp.pi/2); PIO2_LS = f32(mp.pi/2 - mp.mpf(PIO2_HS)); PI_HS = f32(mp.pi)
PIO4_HS = f32(mp.pi/4)
def asin_R32(t):
    p = f32(t * horner32(PS_S, t))
    q = f32(1.0 + f32(t * QS1_S))
    return f32(p / q)
def asin_f32(x):
    x = f32(x)
    if x != x: return float('nan')
    ax = f32(abs(x)); sgn = bits32(x) >> 31
    if ax > 1.0: return float('nan')
    if ax == 1.0: return f32(f32(x * PIO2_HS) + f32(x * PIO2_LS))
    if ax < 0.5:
        if ax < 2.44e-4: return x
        return f32(x + f32(x * asin_R32(f32(x * x))))
    w = f32(1.0 - ax); t = f32(w * 0.5)
    s = f32(math.sqrt(t)); r = asin_R32(t)
    if ax >= 0.975:
        res = f32(PIO2_HS - f32(2.0 * f32(s + f32(s * r))))
    else:
        df = of_bits32(bits32(s) & 0xfffff000)
        c  = f32(f32(t - f32(df * df)) / f32(s + df))
        p2 = f32(f32(2.0 * f32(s * r)) - f32(PIO2_LS - f32(2.0 * c)))
        q2 = f32(PIO4_HS - f32(2.0 * df))
        res = f32(PIO4_HS - f32(p2 - q2))
    return f32(-res) if sgn else res
def acos_f32(x):
    x = f32(x)
    if x != x: return float('nan')
    ax = f32(abs(x)); neg = bits32(x) >> 31
    if ax > 1.0: return float('nan')
    if ax == 1.0: return 0.0 if not neg else f32(PI_HS + f32(2.0 * PIO2_LS))
    if ax < 0.5:
        if ax < 1.49e-8: return PIO2_HS                # |x| < 2^-26
        z = f32(x * x); r = asin_R32(z)
        return f32(PIO2_HS - f32(x - f32(PIO2_LS - f32(x * r))))
    if neg:
        z = f32(f32(1.0 + x) * 0.5); s = f32(math.sqrt(z)); r = asin_R32(z)
        w = f32(f32(r * s) - PIO2_LS)
        return f32(PI_HS - f32(2.0 * f32(s + w)))
    z = f32(f32(1.0 - x) * 0.5); s = f32(math.sqrt(z)); r = asin_R32(z)
    w = f32(f32(r * s) + 0.0)                      # f32: no head/tail split
    # use simple form: acos = 2*asin(sqrt((1-x)/2)) = 2*(s + s*r)
    return f32(2.0 * f32(s + f32(s * r)))

def check_ainv_rs(which):
    fn = asin_f32 if which == 'asin' else acos_f32
    tru = mp.asin if which == 'asin' else mp.acos
    print(f"# {which} @rs: fdlibm rational kernel (f32)")
    if which == 'asin':
        print("PS: " + " ".join(hexs(c) for c in PS_S) + "  QS1: " + hexs(QS1_S))
        print(f"PIO2_H={hexs(PIO2_HS)} PIO2_L={hexs(PIO2_LS)} PI_H={hexs(PI_HS)}")
    worst = 0.0; xw = None
    for t in range(-1000000, 1000001, 7):
        x = f32(mp.mpf(t) / 1000000); got = fn(x); tr = tru(mp.mpf(x))
        if abs(tr) < 1e-7 or not math.isfinite(got): continue
        e = abs(ulps32(got, tr))
        if e > worst: worst, xw = e, x
    print(f"max error over x in [-1,1]: {worst:.3f} ULP  at x={xw}")
    for x in [0.0, 0.5, 1.0, -1.0, 0.25, 0.75, 0.9, -0.6, 0.1]:
        print(f"  {which}({x}) -> {hexs(fn(x))}   in={hexs(x)}")
    for name, x in [('nan', float('nan')), ('1.5', 1.5)]:
        print(f"  {which}({name}) -> {hexs(fn(x))}")

# ---- tan @rs: 3-part pi/2 reduction + f32 __kernel_tan ----
TAN_T_S = [f32(c) for c in TAN_T[:7]]
PIO4_S = f32(mp.pi/4); PIO4LO_S = f32(mp.pi/4 - mp.mpf(PIO4_S))
TAN_BIG_S = f32(0.6744)
def head0_32(x): return of_bits32(bits32(x) & 0xfffff000)
def ktan32(x, y, iy):
    hx_neg = bits32(x) >> 31
    big = f32(abs(x)) >= TAN_BIG_S
    if big:
        if hx_neg: x = f32(-x); y = f32(-y)
        z = f32(PIO4_S - x); w0 = f32(PIO4LO_S - y); x = f32(z + w0); y = 0.0
    z = f32(x*x); w = f32(z*z)
    r = horner32([TAN_T_S[1],TAN_T_S[3],TAN_T_S[5]], w)
    v = f32(z * horner32([TAN_T_S[2],TAN_T_S[4],TAN_T_S[6]], w))
    s = f32(z * x)
    r = f32(y + f32(z * f32(f32(s * f32(r + v)) + y)))
    r = f32(r + f32(TAN_T_S[0] * s))
    w = f32(x + r)
    if big:
        fac = 1.0 if not hx_neg else -1.0
        v = float(iy)
        return f32(fac * f32(v - f32(2.0 * f32(x - f32(f32(f32(w*w) / f32(w+v)) - r)))))
    if iy == 1: return w
    z2 = head0_32(w); v2 = f32(r - f32(z2 - x)); a = f32(-1.0 / w); t2 = head0_32(a)
    s2 = f32(1.0 + f32(t2 * z2))
    return f32(t2 + f32(a * f32(s2 + f32(t2 * v2))))
def tan_f32(x):
    x = f32(x)
    if x != x or x == INF or x == -INF: return float('nan')
    if x == 0.0: return x
    neg = bits32(x) >> 31; ax = f32(abs(x))
    q, rhi, rlo = reduce_pio2_32(ax)
    iy = -1 if (q & 1) else 1
    t = ktan32(rhi, rlo, iy)
    return f32(-t) if neg else t

def old_ratio_tan_f32(x):
    """math.hoon's OLD @rs +tan (superseded 2026-07-03): (div (sin x) (cos
    x)).  Kept only as the historical comparison point (~1.2 ULP) in the
    check output below -- @rs now has its own dedicated kernel, see
    ktan32_shipped."""
    x = f32(x)
    if x != x or x == INF or x == -INF: return float('nan')
    if x == 0.0: return x
    s = sin_f32(x); c = cos_f32(x)
    if c == 0.0: return float('nan')
    return f32(s / c)

_RS_TANQ = [of_bits32(b) for b in
            (0x3b7bdd30, 0x3a9b5ccb, 0x3c233846, 0x3cb11eca,
             0x3d5d294d, 0x3e088844, 0x3eaaaaab)]

def ktan32_shipped(rhi, rlo, iy_is_tan):
    """math.hoon's ACTUAL @rs +tan kernel (shipped 2026-07-03, mirrored in
    the C jet _rs_tan): Q(z)=(tan(r)/r-1)/z, Chebyshev-fit (7 coeffs,
    degree 6) evaluated via Horner, with the dominant linear term rhi added
    LAST (w2=rhi+r) rather than multiplied through the polynomial -- an
    earlier abandoned draft's failure (9.68 ULP, worse than the ratio it
    was meant to replace) traced to getting this ordering wrong, not to
    f32 precision being insufficient.  iy_is_tan=False takes the -cot path
    (odd reduction quadrant) via a plain reciprocal."""
    z = f32(rhi * rhi)
    acc = _RS_TANQ[0]
    for coeff in _RS_TANQ[1:]:
        acc = f32(f32(acc * z) + coeff)
    s = f32(z * rhi)
    corrpoly = f32(s * acc)
    rlo_corr = f32(rlo * f32(1.0 + z))
    r_total = f32(corrpoly + rlo_corr)
    w2 = f32(rhi + r_total)
    if iy_is_tan: return w2
    return f32(-1.0 / w2)

def tan_f32_shipped(x):
    x = f32(x)
    if x != x or x == INF or x == -INF: return float('nan')
    if x == 0.0: return x
    neg = bits32(x) >> 31; ax = f32(abs(x))
    q, rhi, rlo = reduce_pio2_32(ax)
    t = ktan32_shipped(rhi, rlo, (q & 1) == 0)
    return f32(-t) if neg else t

def check_tan_rs():
    print("# tan @rs: dedicated kernel (shipped 2026-07-03; see ktan32_shipped), q*pi/2 reduction + Chebyshev-fit Q(z), -cot path for odd quadrants")
    ws = 0.0; xs = None      # shipped: what math.hoon / the C jet actually compute
    wo = 0.0; xo = None      # old composed sin/cos ratio, historical comparison only
    for t in range(-200000, 200001, 7):
        x = f32(mp.mpf(t) / 1000); tr = mp.tan(mp.mpf(x))
        if not math.isfinite(float(tr)) or abs(tr) > 1e7 or abs(tr) < 1e-6: continue
        gs = tan_f32_shipped(x)
        if math.isfinite(gs):
            e = abs(ulps32(gs, tr))
            if e > ws: ws, xs = e, x
        go = old_ratio_tan_f32(x)
        if math.isfinite(go):
            e = abs(ulps32(go, tr))
            if e > wo: wo, xo = e, x
    print(f"  max {ws:.3f} ULP at x={xs}  (dedicated kernel, shipped; sampled)")
    print(f"  [historical comparison, superseded] old sin/cos ratio: max {wo:.3f} ULP at x={xo}")
    for x in [0.0, 0.5, 1.0, -1.0, 0.7853982, 2.0, 10.0, 100.0]:
        print(f"  tan({x}) -> {hexs(tan_f32_shipped(x))}   in={hexs(x)}")
    for name, x in [('+inf', INF), ('nan', float('nan')), ('-0', -0.0)]:
        o = tan_f32_shipped(x)
        ib = "0x7f800000" if x==INF else "0x7fc00000" if x!=x else hexs(x)
        print(f"  tan({name}) -> {hexs(o)}   in={ib}")

# ==================== @rd / @rs composed ops (sqt, cbrt, pow, pow-n, log-2, log-10, atan2) ====================
#  These arms in /lib/math are themselves literal compositions of +exp/+log/+atan
#  plus a few extra scalar constants (see math.hoon): pow(x,n) = exp(n*log x)
#  (or repeated-multiply pow-n for positive integer n), cbrt(x) = sign(x)*
#  exp(log|x|/3), log-2/log-10 = e*log_b(2) + log(m)/ln(b) via the +lr mantissa/
#  exponent split, atan2 = quadrant dispatch over +atan.  Rather than a
#  separate exact-arithmetic oracle, this section SIMULATES the actual Hoon
#  composition using the already-validated exp_f64/log_f64/atan_f64 (deg-11/
#  deg-9/deg-10 minimax, confirmed bit-identical to math.hoon's baked hex
#  above) plus math.sqrt/math.fma (both correctly-rounded IEEE ops) for +sqt's
#  Markstein step -- so this is a real re-execution of the algorithm, not a
#  derived bound, for every function below EXCEPT where noted.
EXP_COEFFS_RD = gen_exp_coeffs(11)
LOG_COEFFS_RD = gen_log_coeffs(9)
EXP_COEFFS_RS = gen_exp_coeffs_f32(6)
LOG_COEFFS_RS = gen_log_coeffs_f32(4)

def exp64(x):  return exp_f64(x, EXP_COEFFS_RD)
def log64(x):  return log_f64(x, LOG_COEFFS_RD)
def exp32(x):  return exp_f32(x, EXP_COEFFS_RS)
def log32(x):  return log_f32(x, LOG_COEFFS_RS)

def lr_f64(x):
    """[[ef, l1] mirrors math.hoon's +lr:rd -- same reduction as log_f64, but
    returns the un-combined [exponent-as-@rd, log(mantissa)] pair (needed by
    +log-2/+log-10 so the integer part is added with no division rounding)."""
    x = f64(x)
    sub = x < 2.2250738585072014e-308
    xx = f64(x * 18014398509481984.0) if sub else x
    ae = -54 if sub else 0
    b = bits(xx); ef = ((b >> 52) & 0x7ff) - 1023
    m = of_bits((b & ((1<<52)-1)) | 0x3ff0000000000000)
    if m >= SQRT2: m = f64(m * 0.5); ef += 1
    ef += ae
    f = f64(m - ONE)
    s = f64(f / f64(m + ONE))
    z = f64(s * s)
    p2 = horner(LOG_COEFFS_RD, z)
    r = f64(f64(z + z) * p2)
    l1 = f64(f - f64(s * f64(f - r)))
    return f64(float(ef)), l1

def lr_f32(x):
    x = f32(x)
    sub = x < 1.1754943508222875e-38
    xx = f32(x * 16777216.0) if sub else x
    ae = -24 if sub else 0
    b = bits32(xx); ef = ((b >> 23) & 0xff) - 127
    m = struct.unpack('>f', struct.pack('>I', (b & 0x7fffff) | 0x3f800000))[0]
    if m >= SQRT2_S: m = f32(m * 0.5); ef += 1
    ef += ae
    f = f32(m - 1.0)
    s = f32(f / f32(m + 1.0))
    z = f32(s * s)
    p2 = horner32(LOG_COEFFS_RS, z)
    r = f32(f32(z + z) * p2)
    l1 = f32(f - f32(s * f32(f - r)))
    return f32(float(ef)), l1

LOG10_2_RD = f64('0.30102999566398119521373889472449')   # log10(2)
INVLN10_RD = f64('0.43429448190325182765112891891661')   # 1/ln10
LOG10_2_RS = f32('0.30102999566398119521373889472449')
INVLN10_RS = f32('0.43429448190325182765112891891661')

def log2_f64(x):
    ef, lm = lr_f64(x)
    return f64(ef + f64(lm * LOG2E))
def log10_f64(x):
    ef, lm = lr_f64(x)
    return f64(f64(ef * LOG10_2_RD) + f64(lm * INVLN10_RD))
def log2_f32(x):
    ef, lm = lr_f32(x)
    return f32(ef + f32(lm * LOG2E_S))
def log10_f32(x):
    ef, lm = lr_f32(x)
    return f32(f32(ef * LOG10_2_RS) + f32(lm * INVLN10_RS))

def sqt_f64(x):
    """Markstein-corrected f64 sqrt (mirrors +sqt:rd).  math.sqrt is itself
    IEEE-mandated correctly-rounded, so this seed is already exact; the
    Markstein step (which in production corrects a merely-faithful native
    sqrt) is therefore a no-op here -- consistent with the Hoon comment's own
    claim that the composed result is correctly rounded."""
    x = f64(x)
    if x != x: return float('nan')
    if x == INF: return INF
    if x == 0.0: return x
    if x < 0.0: return float('nan')
    g = f64(math.sqrt(x))
    h = f64(0.5 / g)
    r = f64(math.fma(-g, g, x))
    return f64(math.fma(h, r, g))
def sqt_f32(x):
    """+sqt:rs delegates directly to the native (correctly-rounded) f32 sqrt;
    modeled as double-precision sqrt rounded once to f32 (safe: 53 vs 24
    mantissa bits is far more headroom than the well-known f64->f32
    double-rounding-safety margin)."""
    x = f32(x)
    if x != x: return float('nan')
    if x == INF: return INF
    if x == 0.0: return x
    if x < 0.0: return float('nan')
    return f32(math.sqrt(x))

THIRD_RD = f64(mp.mpf(1) / 3)
THIRD_RS = f32(mp.mpf(1) / 3)
def cbt_f64(x):
    x = f64(x)
    if x != x: return x
    if x == 0.0: return x
    ax = f64(abs(x))
    r = exp64(f64(log64(ax) * THIRD_RD))
    return f64(-r) if bits(x) >> 63 else r
def cbt_f32(x):
    x = f32(x)
    if x != x: return x
    if x == 0.0: return x
    ax = f32(abs(x))
    r = exp32(f32(log32(ax) * THIRD_RS))
    return f32(-r) if bits32(x) >> 31 else r

def pow_n_f64(x, n):
    n = int(n)
    if n == 0: return 1.0
    p = x
    for _ in range(n - 1): p = f64(p * x)
    return p
def pow_n_f32(x, n):
    n = int(n)
    if n == 0: return 1.0
    p = x
    for _ in range(n - 1): p = f32(p * x)
    return p
def pow_f64(x, n):
    if n == f64(round(n)) and n > 0:
        return pow_n_f64(x, int(n))
    return exp64(f64(n * log64(x)))
def pow_f32(x, n):
    if n == f32(round(n)) and n > 0:
        return pow_n_f32(x, int(n))
    return exp32(f32(n * log32(x)))

PI_RD = f64(mp.pi); PI_RS = f32(mp.pi)
def atan2_f64(y, x):
    y, x = f64(y), f64(x)
    if x > 0.0: return atan_f64(f64(y / x))
    if x < 0.0 and y >= 0.0: return f64(atan_f64(f64(y / x)) + PI_RD)
    if x < 0.0 and y < 0.0: return f64(atan_f64(f64(y / x)) - PI_RD)
    if x == 0.0 and y > 0.0: return f64(PI_RD / 2)
    if x == 0.0 and y < 0.0: return f64(-1.0 * f64(PI_RD / 2))
    return 0.0
def atan2_f32(y, x):
    y, x = f32(y), f32(x)
    if x > 0.0: return atan_f32(f32(y / x))
    if x < 0.0 and y >= 0.0: return f32(atan_f32(f32(y / x)) + PI_RS)
    if x < 0.0 and y < 0.0: return f32(atan_f32(f32(y / x)) - PI_RS)
    if x == 0.0 and y > 0.0: return f32(PI_RS / 2)
    if x == 0.0 and y < 0.0: return f32(-1.0 * f32(PI_RS / 2))
    return 0.0

def _sweep_1d(fn, truefn, xs, skip_small=0.0, ulpfn=None, roundfn=None):
    worst = 0.0; xw = None; n = 0
    for xv in xs:
        x = roundfn(xv) if roundfn else xv
        got = fn(x)
        if not math.isfinite(got): continue
        tru = truefn(mp.mpf(x) if not isinstance(x, mp.mpf) else x)
        if abs(tru) < skip_small: continue
        e = abs(ulpfn(got, tru))
        n += 1
        if e > worst: worst, xw = e, x
    return worst, xw, n

def check_log2log10_rd():
    print("# log-2/log-10 @rd: +lr mantissa/exponent split + linear combine (composed; full simulation)")
    for name, fn, truefn in [('log-2', log2_f64, lambda xm: mp.log(xm)/mp.log(2)),
                              ('log-10', log10_f64, lambda xm: mp.log(xm)/mp.log(10))]:
        xs = [f64(mp.mpf(t)/1000) for t in range(1, 200001, 7)]
        worst, xw, n = _sweep_1d(fn, truefn, xs, ulpfn=ulps)
        print(f"  {name}: max {worst:.3f} ULP at x={xw}  (sampled, n={n}, x in (0,200])")
        for v in [1.0, 2.0, 0.5, 10.0, 100.0, 0.1]:
            print(f"    {name}({v}) -> {hexd(fn(v))}")

def check_log2log10_rs():
    print("# log-2/log-10 @rs: +lr mantissa/exponent split + linear combine (composed; full simulation)")
    for name, fn, truefn in [('log-2', log2_f32, lambda xm: mp.log(xm)/mp.log(2)),
                              ('log-10', log10_f32, lambda xm: mp.log(xm)/mp.log(10))]:
        xs = [f32(mp.mpf(t)/1000) for t in range(1, 200001, 7)]
        worst, xw, n = _sweep_1d(fn, truefn, xs, ulpfn=ulps32)
        print(f"  {name}: max {worst:.3f} ULP at x={xw}  (sampled, n={n}, x in (0,200])")

def check_sqt_rd():
    print("# sqt @rd: Markstein-corrected sqrt (composed; full simulation)")
    xs = [f64(mp.mpf(t)/1000) for t in range(1, 2000001, 37)]
    worst, xw, n = _sweep_1d(sqt_f64, mp.sqrt, xs, ulpfn=ulps)
    print(f"  max {worst:.3f} ULP at x={xw}  (sampled, n={n}, x in (0,2000])")
def check_sqt_rs():
    print("# sqt @rs: native f32 sqrt (composed; full simulation)")
    xs = [f32(mp.mpf(t)/1000) for t in range(1, 2000001, 37)]
    worst, xw, n = _sweep_1d(sqt_f32, mp.sqrt, xs, ulpfn=ulps32)
    print(f"  max {worst:.3f} ULP at x={xw}  (sampled, n={n}, x in (0,2000])")

def _cbrt_true(xm):                    # real cube root, defined for negative reals too
    return -mp.cbrt(-xm) if xm < 0 else mp.cbrt(xm)

def check_cbrt_rd():
    print("# cbrt @rd: sign(x)*exp(log|x|/3) (composed; full simulation)")
    xs = [f64(mp.mpf(t)/1000) for t in range(-2000000, 2000001, 37) if t != 0]
    worst, xw, n = _sweep_1d(cbt_f64, _cbrt_true, xs, ulpfn=ulps)
    print(f"  max {worst:.3f} ULP at x={xw}  (sampled, n={n}, x in [-2000,2000])")
def check_cbrt_rs():
    print("# cbrt @rs: sign(x)*exp(log|x|/3) (composed; full simulation)")
    xs = [f32(mp.mpf(t)/1000) for t in range(-2000000, 2000001, 37) if t != 0]
    worst, xw, n = _sweep_1d(cbt_f32, _cbrt_true, xs, ulpfn=ulps32)
    print(f"  max {worst:.3f} ULP at x={xw}  (sampled, n={n}, x in [-2000,2000])")

def check_pow_rd():
    print("# pow/pow-n @rd: pow-n repeated-multiply (int n>0) else exp(n*log x) (composed; full simulation)")
    worst = 0.0; xw = None; n = 0
    for xt in range(1, 2001, 3):
        x = f64(mp.mpf(xt)/100)
        for nt in range(-2000, 2001, 11):
            nn = f64(mp.mpf(nt)/100)
            got = pow_f64(x, nn)
            if not math.isfinite(got): continue
            tru = mp.mpf(x) ** mp.mpf(nn)
            if tru == 0 or not math.isfinite(float(tru)): continue
            e = abs(ulps(got, tru)); n += 1
            if e > worst: worst, xw = e, (x, nn)
    print(f"  max {worst:.3f} ULP at (x,n)={xw}  (sampled grid, n={n}, x in (0,20], n in [-20,20])")
    print("  pow-n exact-int-exponent spot check (n=2,3,5,-3):")
    for nn in [2, 3, 5, -3]:
        worst2 = 0.0; xw2 = None; m = 0
        for xt in range(1, 3001, 7):
            x = f64(mp.mpf(xt)/100)
            got = pow_f64(x, float(nn)) if nn > 0 else exp64(f64(float(nn)*log64(x)))
            tru = mp.mpf(x) ** nn
            if not math.isfinite(got): continue
            e = abs(ulps(got, tru)); m += 1
            if e > worst2: worst2, xw2 = e, x
        print(f"    n={nn}: max {worst2:.3f} ULP at x={xw2} (n={m})")

def check_pow_rs():
    print("# pow/pow-n @rs: composed; full simulation")
    worst = 0.0; xw = None; n = 0
    for xt in range(1, 2001, 3):
        x = f32(mp.mpf(xt)/100)
        for nt in range(-2000, 2001, 11):
            nn = f32(mp.mpf(nt)/100)
            got = pow_f32(x, nn)
            if not math.isfinite(got): continue
            tru = mp.mpf(x) ** mp.mpf(nn)
            if tru == 0 or not math.isfinite(float(tru)): continue
            e = abs(ulps32(got, tru)); n += 1
            if e > worst: worst, xw = e, (x, nn)
    print(f"  max {worst:.3f} ULP at (x,n)={xw}  (sampled grid, n={n})")

def check_atan2_rd():
    print("# atan2 @rd: quadrant dispatch over +atan (composed; full simulation)")
    worst = 0.0; xw = None; n = 0
    for yt in range(-2000, 2001, 13):
        for xt in range(-2000, 2001, 13):
            y = f64(mp.mpf(yt)/100); x = f64(mp.mpf(xt)/100)
            if y == 0.0 and x == 0.0: continue
            got = atan2_f64(y, x)
            tru = mp.atan2(mp.mpf(y), mp.mpf(x))
            if not math.isfinite(got): continue
            e = abs(ulps(got, tru)); n += 1
            if e > worst: worst, xw = e, (y, x)
    print(f"  max {worst:.3f} ULP at (y,x)={xw}  (sampled grid, n={n}, [-20,20]x[-20,20])")
def check_atan2_rs():
    print("# atan2 @rs: composed; full simulation")
    worst = 0.0; xw = None; n = 0
    for yt in range(-2000, 2001, 13):
        for xt in range(-2000, 2001, 13):
            y = f32(mp.mpf(yt)/100); x = f32(mp.mpf(xt)/100)
            if y == 0.0 and x == 0.0: continue
            got = atan2_f32(y, x)
            tru = mp.atan2(mp.mpf(y), mp.mpf(x))
            if not math.isfinite(got): continue
            e = abs(ulps32(got, tru)); n += 1
            if e > worst: worst, xw = e, (y, x)
    print(f"  max {worst:.3f} ULP at (y,x)={xw}  (sampled grid, n={n})")

# ============================ @rh (f16) ============================
#  Unlike @rd/@rs above (which independently DERIVE minimax coefficients via
#  mp.chebyfit and were confirmed bit-identical to math.hoon's baked hex --
#  see the exp@rd cross-check), @rh's polynomials below are DECODED DIRECTLY
#  from the literal hex constants in libmath/desk/lib/math.hoon (the +exp/
#  +log/+lr/+sin/+cos/+atan/+asin/+acos arms under `++  rh`) rather than
#  regenerated, since @rh's degree/breakpoints were hand-tuned and an
#  independent re-fit is not guaranteed to reproduce the exact same bits.
#  This means the @rh checks below test the ACTUAL shipped algorithm bit for
#  bit, which is the strongest possible oracle.  Exhaustive: iterates all
#  2^16 @rh bit patterns, per the paper's own claim for half precision.
def f16(x):
    x = float(x)
    try:
        return struct.unpack('>e', struct.pack('>e', x))[0]
    except (OverflowError, struct.error):
        # struct's 'e' packer raises instead of rounding-to-infinity on
        # overflow (unlike hardware IEEE conversion); match IEEE round-to-
        # nearest-even-ties overflow-to-infinity semantics by hand.  f16 max
        # finite is 65504.0; halfway to the next (nonexistent) representable
        # value is 65520.0, so ties-to-even rounds 65504.0..65520.0 down.
        if x != x: return float('nan')
        if abs(x) < 65520.0: return math.copysign(65504.0, x)
        return math.copysign(math.inf, x)
def bits16(x):return struct.unpack('>H', struct.pack('>e', f16(x)))[0]
def hex16(x): return f"0x{bits16(x):04x}"
def of_bits16(b): return struct.unpack('>e', struct.pack('>H', b))[0]

def ulps16(approx, true_mpf):
    """Local ULP spacing computed via bit-pattern increment (NOT
    math.nextafter on the double representation followed by re-rounding to
    f16 -- that double-rounds to a no-op near/below the smallest subnormal,
    since a double-precision nudge is far too small to survive a second
    round-to-f16, which silently produced a bogus ulp of a*2^-10 there;
    caught by exp@rh showing a spurious ~500 ULP 'error' at a value that
    was actually already correctly rounded to the nearest subnormal)."""
    a = f16(approx)
    if a != a or not math.isfinite(a): return float('inf')
    ab = bits16(abs(a)) & 0x7fff
    up_mag = of_bits16(ab + 1) if ab + 1 < 0x7c00 else math.inf
    ulp = up_mag - abs(a)
    if ulp == 0: ulp = 2 ** -24
    return float((mp.mpf(a) - true_mpf) / mp.mpf(ulp))

def horner16(coeffs, r):              # ascending coeffs; strict-f16 Horner (== Hoon's roll-of-flop)
    acc = f16(coeffs[-1])
    for c in reversed(coeffs[:-1]):
        acc = f16(f16(acc * r) + c)
    return acc

def _sweep_rh_exhaustive(fn, truefn, skip_small=0.0, skip_nan_out=True):
    worst = 0.0; xw = None; n = 0
    for b in range(65536):
        x = of_bits16(b)
        if x != x: continue
        got = fn(x)
        if skip_nan_out and (got != got or not math.isfinite(got)): continue
        tru = truefn(mp.mpf(x))
        if not mp.isfinite(tru) or abs(tru) < skip_small: continue
        e = abs(ulps16(got, tru))
        n += 1
        if e > worst: worst, xw = e, x
    return worst, xw, n

# ---- exp @rh (Cody-Waite + degree-4 minimax; math.hoon ++exp:rh) ----
LOG2E_H = of_bits16(0x3dc5); LN2HI_H = of_bits16(0x3980); LN2LO_H = of_bits16(0x1dc8)
EXP_C_H = [of_bits16(h) for h in (0x3c00, 0x3c00, 0x3800, 0x3160, 0x295c)]
def exp_f16(x):
    x = f16(x)
    if x != x: return float('nan')
    if bits16(x) == 0x7c00: return math.inf
    if bits16(x) == 0xfc00: return 0.0
    xk = f16(x * LOG2E_H)
    if not math.isfinite(xk): return math.inf if xk > 0 else 0.0
    k = int(math.floor(xk + 0.5))
    if k >= 17: return math.inf
    if k <= -25: return 0.0
    kf = f16(float(k))
    r = f16(f16(x - f16(kf * LN2HI_H)) - f16(kf * LN2LO_H))
    p = horner16(EXP_C_H, r)
    try:    return f16(math.ldexp(p, k))
    except OverflowError: return math.inf

def check_exp_rh():
    print("# exp @rh: Cody-Waite reduction + degree-4 minimax poly (exact math.hoon hex; exhaustive)")
    worst, xw, n = _sweep_rh_exhaustive(exp_f16, mp.exp)
    print(f"  max {worst:.3f} ULP at x={xw}  (EXHAUSTIVE over 2^16, n={n})")

# ---- log/lr @rh (x=2^e*m reduction + degree-1 atanh; math.hoon ++log/++lr:rh) ----
LOG_C_H = [of_bits16(h) for h in (0x3555, 0x3266)]
def lr_f16(x):                        # x: finite positive @rh -> (ef, l1)
    b = bits16(x)
    sub = ((b >> 10) & 0x1f) == 0
    xx = f16(x * of_bits16(0x6400)) if sub else x
    ae = -10 if sub else 0
    b2 = bits16(xx)
    e = ((b2 >> 10) & 0x1f) - 15
    m = of_bits16((b2 & 0x3ff) | 0x3c00)
    big = m >= of_bits16(0x3da8)
    if big: m = f16(m * of_bits16(0x3800)); e += 1
    e += ae
    f = f16(m - of_bits16(0x3c00))
    s = f16(f / f16(m + of_bits16(0x3c00)))
    z = f16(s * s)
    p2 = horner16(LOG_C_H, z)
    r = f16(f16(z + z) * p2)
    l1 = f16(f - f16(s * f16(f - r)))
    return f16(float(e)), l1

def log_f16(x):
    x = f16(x)
    if x != x: return float('nan')
    if bits16(x) == 0x7c00: return math.inf
    if bits16(x) in (0x0, 0x8000): return -math.inf
    if bits16(x) >> 15 == 1: return float('nan')
    ef, l1 = lr_f16(x)
    hi = f16(ef * LN2HI_H); lo = f16(ef * LN2LO_H)
    return f16(hi + f16(l1 + lo))

LOG10_2_H = of_bits16(0x34d1); INVLN10_H = of_bits16(0x36f3)
def log2_f16(x):
    x = f16(x)
    if x != x: return float('nan')
    if bits16(x) == 0x7c00: return math.inf
    if bits16(x) in (0x0, 0x8000): return -math.inf
    if bits16(x) >> 15 == 1: return float('nan')
    ef, lm = lr_f16(x)
    return f16(ef + f16(lm * LOG2E_H))
def log10_f16(x):
    x = f16(x)
    if x != x: return float('nan')
    if bits16(x) == 0x7c00: return math.inf
    if bits16(x) in (0x0, 0x8000): return -math.inf
    if bits16(x) >> 15 == 1: return float('nan')
    ef, lm = lr_f16(x)
    return f16(f16(ef * LOG10_2_H) + f16(lm * INVLN10_H))

def check_log_rh():
    print("# log @rh: x=2^e*m reduction + degree-1 atanh poly (exact math.hoon hex; exhaustive)")
    worst, xw, n = _sweep_rh_exhaustive(log_f16, mp.log)
    print(f"  max {worst:.3f} ULP at x={xw}  (EXHAUSTIVE over 2^16 (x>0), n={n})")
def check_log2log10_rh():
    print("# log-2/log-10 @rh: +lr split + linear combine (exact math.hoon hex; exhaustive)")
    for name, fn, truefn in [('log-2', log2_f16, lambda xm: mp.log(xm)/mp.log(2)),
                              ('log-10', log10_f16, lambda xm: mp.log(xm)/mp.log(10))]:
        worst, xw, n = _sweep_rh_exhaustive(fn, truefn)
        print(f"  {name}: max {worst:.3f} ULP at x={xw}  (EXHAUSTIVE over 2^16 (x>0), n={n})")

# ---- sin/cos/tan @rh (quarter-turn reduction + 2-coeff kernels; math.hoon ++rh-trig) ----
SC_H = [of_bits16(h) for h in (0xb155, 0x2044)]
CC_H = [of_bits16(h) for h in (0x2955, 0x95b0)]
def ksin16(xx, yy):
    z = f16(xx * xx)
    r = horner16(SC_H[1:], z)
    v = f16(z * xx)
    aa = f16(f16(of_bits16(0x3800) * yy) - f16(v * r))
    bb = f16(f16(z * aa) - yy)
    dd = f16(bb - f16(v * SC_H[0]))
    return f16(xx - dd)
def kcos16(xx, yy):
    z = f16(xx * xx)
    rc = horner16(CC_H, z)
    hz = f16(of_bits16(0x3800) * z)
    w2 = f16(of_bits16(0x3c00) - hz)
    aa = f16(f16(of_bits16(0x3c00) - w2) - hz)
    bb = f16(f16(f16(z * z) * rc) - f16(xx * yy))
    return f16(w2 + f16(aa + bb))
def trig_fin16(is_sin, ax, sb):
    q = int(math.floor(f16(ax * of_bits16(0x3918)) + 0.5))
    qf = f16(float(abs(q)))
    r1 = f16(ax - f16(qf * of_bits16(0x3e00)))
    r2 = f16(r1 - f16(qf * of_bits16(0x2c80)))
    w = f16(qf * of_bits16(0xfed))
    rhi = f16(r2 - w); rlo = f16(f16(r2 - rhi) - w)
    m = abs(q) & 3
    ks = ksin16(rhi, rlo); kc = kcos16(rhi, rlo)
    if is_sin:
        v = ks if m == 0 else kc if m == 1 else (f16(-ks) if m == 2 else f16(-kc))
        return f16(-v) if sb == 1 else v
    return kc if m == 0 else (f16(-ks) if m == 1 else (f16(-kc) if m == 2 else ks))
def sin_f16(x):
    x = f16(x)
    if x != x: return float('nan')
    if bits16(x) in (0x7c00, 0xfc00): return float('nan')
    if bits16(x) in (0x0, 0x8000): return x
    return trig_fin16(True, of_bits16(bits16(x) & 0x7fff), bits16(x) >> 15)
def cos_f16(x):
    x = f16(x)
    if x != x: return float('nan')
    if bits16(x) in (0x7c00, 0xfc00): return float('nan')
    return trig_fin16(False, of_bits16(bits16(x) & 0x7fff), 0)
def tan_f16(x): return f16(sin_f16(x) / cos_f16(x))

def check_trig_rh(which):
    fn = sin_f16 if which == 'sin' else cos_f16
    tru = mp.sin if which == 'sin' else mp.cos
    print(f"# {which} @rh: quarter-turn reduction + 2-coeff kernel (exact math.hoon hex; exhaustive)")
    worst, xw, n = _sweep_rh_exhaustive(fn, tru, skip_small=1e-3)
    print(f"  max {worst:.3f} ULP at x={xw}  (EXHAUSTIVE over 2^16 (|f|>1e-3), n={n})")
def check_tan_rh():
    print("# tan @rh: sin/cos ratio (exact math.hoon hex; exhaustive)")
    worst, xw, n = _sweep_rh_exhaustive(tan_f16, mp.tan, skip_small=1e-3)
    print(f"  max {worst:.3f} ULP at x={xw}  (EXHAUSTIVE over 2^16 (|f|>1e-3), n={n})")

# ---- atan @rh (fdlibm breakpoint reduction + degree-2 minimax; math.hoon ++atan/++rh-atan) ----
AT_H = [of_bits16(h) for h in (0x3555, 0xb266, 0x3092)]
def atan_ker16(ax):
    one = of_bits16(0x3c00); two = of_bits16(0x4000); ohf = of_bits16(0x3e00)
    if ax < of_bits16(0x3700):
        xr, hi, lo, direct = ax, of_bits16(0), of_bits16(0), True
    elif ax < of_bits16(0x3980):
        xr = f16(f16(f16(ax + ax) - one) / f16(two + ax)); hi = of_bits16(0x376b); lo = of_bits16(0x19c); direct = False
    elif ax < of_bits16(0x3cc0):
        xr = f16(f16(ax - one) / f16(ax + one)); hi = of_bits16(0x3a48); lo = of_bits16(0xbed); direct = False
    elif ax < of_bits16(0x40e0):
        xr = f16(f16(ax - ohf) / f16(one + f16(ohf * ax))); hi = of_bits16(0x3bdd); lo = of_bits16(0x87a1); direct = False
    else:
        xr = f16(of_bits16(0xbc00) / ax); hi = of_bits16(0x3e48); lo = of_bits16(0xfed); direct = False
    z = f16(xr * xr)
    s = f16(z * horner16(AT_H, z))
    if direct: return f16(xr - f16(xr * s))
    return f16(hi - f16(f16(f16(xr * s) - lo) - xr))
def atan_f16(x):
    x = f16(x)
    if x != x: return float('nan')
    if bits16(x) == 0x7c00: return of_bits16(0x3e48)
    if bits16(x) == 0xfc00: return of_bits16(0xbe48)
    if bits16(x) in (0x0, 0x8000): return x
    neg = bits16(x) >> 15
    r = atan_ker16(of_bits16(bits16(x) & 0x7fff))
    return f16(-r) if neg else r
def check_atan_rh():
    print("# atan @rh: fdlibm breakpoint reduction + degree-2 minimax (exact math.hoon hex; exhaustive)")
    worst, xw, n = _sweep_rh_exhaustive(atan_f16, mp.atan, skip_small=1e-4)
    print(f"  max {worst:.3f} ULP at x={xw}  (EXHAUSTIVE over 2^16, n={n})")

# ---- asin/acos @rh (fdlibm rational kernel; math.hoon ++asin/++acos/++rh-ainv) ----
PS_H = [of_bits16(h) for h in (0x3155, 0x2cea, 0x2729, 0x2ccc)]
def rr16(t): return f16(t * horner16(PS_H, t))
def asin_f16(x):
    x = f16(x)
    if x != x: return float('nan')
    sgn = bits16(x) >> 15
    ax = of_bits16(bits16(x) & 0x7fff)
    one = of_bits16(0x3c00)
    if ax > one: return float('nan')
    if ax == one: return f16(f16(x * of_bits16(0x3e48)) + f16(x * of_bits16(0xfed)))
    if ax < of_bits16(0x3800):
        if ax < of_bits16(0xc00): return x
        t = f16(x * x)
        return f16(x + f16(x * rr16(t)))
    w = f16(one - ax); t = f16(w * of_bits16(0x3800))
    r = rr16(t); s = f16(math.sqrt(t))
    if ax >= of_bits16(0x3bcd):
        res = f16(of_bits16(0x3e48) - f16(of_bits16(0x4000) * f16(s + f16(s * r))))
        return f16(-res) if sgn == 1 else res
    df = of_bits16(bits16(s) & 0xfff0)
    c = f16(f16(t - f16(df * df)) / f16(s + df))
    p2 = f16(f16(of_bits16(0x4000) * f16(s * r)) - f16(of_bits16(0xfed) - f16(of_bits16(0x4000) * c)))
    q2 = f16(of_bits16(0x3a48) - f16(of_bits16(0x4000) * df))
    res = f16(of_bits16(0x3a48) - f16(p2 - q2))
    return f16(-res) if sgn == 1 else res
def acos_f16(x):
    x = f16(x)
    if x != x: return float('nan')
    neg = bits16(x) >> 15
    ax = of_bits16(bits16(x) & 0x7fff)
    one = of_bits16(0x3c00)
    if ax > one: return float('nan')
    if ax == one:
        return 0.0 if neg == 0 else f16(of_bits16(0x4248) + f16(of_bits16(0x4000) * of_bits16(0xfed)))
    if ax < of_bits16(0x3800):
        z = f16(x * x); r = rr16(z)
        return f16(of_bits16(0x3e48) - f16(x - f16(of_bits16(0xfed) - f16(x * r))))
    if neg == 1:
        z = f16(f16(one + x) * of_bits16(0x3800)); s = f16(math.sqrt(z)); r = rr16(z)
        w = f16(f16(r * s) - of_bits16(0xfed))
        return f16(of_bits16(0x4248) - f16(of_bits16(0x4000) * f16(s + w)))
    z = f16(f16(one - x) * of_bits16(0x3800)); s = f16(math.sqrt(z)); r = rr16(z)
    return f16(of_bits16(0x4000) * f16(s + f16(s * r)))
def check_ainv_rh(which):
    fn = asin_f16 if which == 'asin' else acos_f16
    tru = mp.asin if which == 'asin' else mp.acos
    print(f"# {which} @rh: fdlibm rational kernel (exact math.hoon hex; exhaustive)")
    def truefn(xm):
        if abs(xm) > 1: return mp.mpf('nan')
        return tru(xm)
    worst, xw, n = _sweep_rh_exhaustive(fn, truefn, skip_small=1e-4)
    print(f"  max {worst:.3f} ULP at x={xw}  (EXHAUSTIVE over 2^16 (|x|<=1), n={n})")

# ---- composed @rh: sqt/cbrt/pow/pow-n/atan2 ----
def sqt_f16(x):
    x = f16(x)
    if x != x: return float('nan')
    if bits16(x) == 0x7c00: return math.inf
    if bits16(x) in (0x0, 0x8000): return x
    if bits16(x) >> 15 == 1: return float('nan')
    return f16(math.sqrt(x))
THIRD_H = of_bits16(0x3555)
def cbt_f16(x):
    x = f16(x)
    if x != x: return x
    if bits16(x) in (0x0, 0x8000): return x
    ax = of_bits16(bits16(x) & 0x7fff)
    r = exp_f16(f16(log_f16(ax) * THIRD_H))
    return f16(-r) if bits16(x) >> 15 else r
def pow_n_f16(x, n):
    n = int(n)
    if n == 0: return 1.0
    p = x
    for _ in range(n - 1): p = f16(p * x)
    return p
def pow_f16(x, n):
    if n == f16(round(n)) and n > 0:
        return pow_n_f16(x, int(n))
    return exp_f16(f16(n * log_f16(x)))
PI_RH = of_bits16(0x4248)              # NB: named PI_RH, not PI_H -- PI_H already denotes @rd's f64 pi (see acos_f64)
def atan2_f16(y, x):
    y, x = f16(y), f16(x)
    if x > 0.0: return atan_f16(f16(y / x))
    if x < 0.0 and y >= 0.0: return f16(atan_f16(f16(y / x)) + PI_RH)
    if x < 0.0 and y < 0.0: return f16(atan_f16(f16(y / x)) - PI_RH)
    if x == 0.0 and y > 0.0: return f16(PI_RH / 2)
    if x == 0.0 and y < 0.0: return f16(-1.0 * f16(PI_RH / 2))
    return 0.0

def check_sqt_rh():
    print("# sqt @rh: native f16 sqrt (exhaustive)")
    worst, xw, n = _sweep_rh_exhaustive(sqt_f16, lambda xm: mp.sqrt(xm) if xm >= 0 else mp.mpf('nan'))
    print(f"  max {worst:.3f} ULP at x={xw}  (EXHAUSTIVE over 2^16 (x>=0), n={n})")
def check_cbrt_rh():
    print("# cbrt @rh: sign(x)*exp(log|x|/3) (exhaustive)")
    worst, xw, n = _sweep_rh_exhaustive(cbt_f16, _cbrt_true)
    print(f"  max {worst:.3f} ULP at x={xw}  (EXHAUSTIVE over 2^16, n={n})")
def check_pow_rh():
    #  NOTE: a naive nested sweep over 2 independent bit-pattern ranges of
    #  2^16 each is O(2^32) -- this is what silently turned an intended
    #  few-second check into an 18+ CPU-minute hang during development.
    #  A ~600x~600 grid (360K pairs, same order as the @rd/@rs pow sweeps
    #  above) is plenty for a composed (non-exhaustive) check.
    print("# pow/pow-n @rh: composed (sampled grid over x,n, NOT exhaustive)")
    worst = 0.0; xw = None; n = 0
    xs = [of_bits16(b) for b in range(1, 31744, 79)]     # positive finite range, ~400 points
    ns = [of_bits16(b) for b in range(0, 65536, 164)]    # ~400 points, full range (+/-)
    for x in xs:
        if x != x or x <= 0 or not math.isfinite(x): continue
        for nn in ns:
            if nn != nn or not math.isfinite(nn): continue
            got = pow_f16(x, nn)
            if not math.isfinite(got): continue
            tru = mp.mpf(x) ** mp.mpf(nn)
            if not mp.isfinite(tru) or tru == 0: continue
            e = abs(ulps16(got, tru)); n += 1
            if e > worst: worst, xw = e, (x, nn)
    print(f"  max {worst:.3f} ULP at (x,n)={xw}  (sampled grid, n={n})")
def check_atan2_rh():
    #  Same O(2^32) trap as check_pow_rh above -- bounded to a ~600x600 grid.
    print("# atan2 @rh: quadrant dispatch over +atan (sampled grid, NOT exhaustive)")
    worst = 0.0; xw = None; n = 0
    ys = [of_bits16(b) for b in range(0, 65536, 164)]    # ~400 points, full range (+/-)
    xs = [of_bits16(b) for b in range(0, 65536, 164)]    # ~400 points, full range (+/-)
    for y in ys:
        if y != y or not math.isfinite(y): continue
        for x in xs:
            if x != x or not math.isfinite(x): continue
            if x == 0.0 and y == 0.0: continue
            got = atan2_f16(y, x)
            if not math.isfinite(got): continue
            tru = mp.atan2(mp.mpf(y), mp.mpf(x))
            e = abs(ulps16(got, tru)); n += 1
            if e > worst: worst, xw = e, (y, x)
    print(f"  max {worst:.3f} ULP at (y,x)={xw}  (sampled grid, n={n})")

# ============================ @rq (f128) ============================
#  numpy float128 is x87 80-bit ext, NOT true IEEE quad, so we model f128 with
#  mpmath rounded to 113-bit RNE per op (= SoftFloat f128 for + - * / ldexp).
#  exp @rq vertical slice; the rest of the rq surface follows the @rd pattern
#  at higher degree.  Run: python3 cheb_check.py exp-rq
QP = 113
def qr(x):                                     # round to f128 (113-bit RNE)
    x = mp.mpf(x)
    if x == 0: return mp.mpf(0)
    m, e = mp.frexp(x); M = mp.nint(m * mp.mpf(2)**QP)
    return M * mp.mpf(2)**(e - QP)
def qhex(v):                                   # 128-bit IEEE-quad encoding
    if v == 0: return "0x" + "0"*32
    s = 1 if v < 0 else 0; a = abs(mp.mpf(v))
    m, e = mp.frexp(a); biased = (e - 1) + 16383
    mant = int(mp.nint((m*2 - 1) * mp.mpf(2)**112))
    if mant == (1 << 112): mant = 0; biased += 1
    return f"0x{(s<<127)|(biased<<112)|mant:032x}"
def qadd(a, b): return qr(mp.mpf(a) + mp.mpf(b))
def qmul(a, b): return qr(mp.mpf(a) * mp.mpf(b))
def qsub(a, b): return qr(mp.mpf(a) - mp.mpf(b))
def qhorner(co, r):
    acc = co[-1]
    for c in reversed(co[:-1]): acc = qadd(qmul(acc, r), c)
    return acc
def qdiv(a, b): return qr(mp.mpf(a) / mp.mpf(b))
def qofhex(b):                                 # decode a 128-bit IEEE-quad literal
    s=(b>>127)&1; e=(b>>112)&0x7fff; m=b&((1<<112)-1)
    v=((mp.mpf(1)+mp.mpf(m)/mp.mpf(2)**112)*mp.mpf(2)**(e-16383)) if e else mp.mpf(m)*mp.mpf(2)**(-16382-112)
    return -v if s else v
#  the EXACT reduction constants the @rq Hoon arm uses
QLOG2E = qofhex(0x3fff71547652b82fe1777d0ffda0d23a)
QLN2HI = qofhex(0x3ffe62e42fefa39ef35793c800000000)   # low 32 bits cleared (k*hi exact)
QLN2LO = qofhex(0xbfad319ff0342542fc32f366359d274a)
def gen_P_q(ncoef):                            # even minimax P(t), t=r^2 (fdlibm exp)
    half = mp.log(2) / 2
    def Pfun(t):
        r = mp.sqrt(t)
        if r == 0: return mp.mpf(1)/6
        E = mp.e**r - 1 - r; c = 2*E/(r+E); return (r-c)/(r*r)
    return [qr(c) for c in reversed(mp.chebyfit(Pfun, [mp.mpf('1e-40'), half*half], ncoef))]
def exp_q(x, P):                               # fdlibm rational reconstruction
    x = qr(x); k = int(mp.floor(qr(x * QLOG2E) + 0.5))
    hi = qsub(x, qmul(mp.mpf(k), QLN2HI)); lo = qmul(mp.mpf(k), QLN2LO); r = qsub(hi, lo)
    t = qmul(r, r); c = qsub(r, qmul(t, qhorner(P, t)))          # c = r - t*P(t)
    y = qsub(mp.mpf(1), qsub(qsub(lo, qdiv(qmul(r, c), qsub(mp.mpf(2), c))), hi))
    return qr(y * mp.mpf(2)**k)
def check_exp_rq():
    P = gen_P_q(11)                            # deg-10 even poly = 11 coeffs
    print("# exp @rq: Cody-Waite + fdlibm rational reconstruction")
    print("#   exp(r) = 1 - ((lo - r*c/(2-c)) - hi),  c = r - t*P(t),  t = r*r")
    print(f"LOG2E={qhex(QLOG2E)}\nLN2HI={qhex(QLN2HI)}\nLN2LO={qhex(QLN2LO)}")
    print("P coeffs (ascending in t=r^2):")
    for i, c in enumerate(P): print(f"  p{i:<2}={qhex(c)}")
    worst = 0
    for t in range(-20000, 20001, 3):          # [-20,20] step 0.003 (rq_check.c+MPFR is authoritative)
        x = qr(mp.mpf(t)/1000); g = exp_q(x, P); tr = mp.e**x
        ulp = mp.mpf(2)**(mp.frexp(g)[1] - QP); worst = max(worst, abs((g-tr)/ulp))
    print(f"max error over [-20,20]: {float(worst):.3f} ULP  (faithful; fdlibm beats the old flat Horner)")
    for v in ['1','0.5','-2','10']:
        print(f"  exp({v}) -> {qhex(exp_q(mp.mpf(v), P))}")

EXP_P_Q = gen_P_q(11)                  # shared by cbrt_q/pow_q (reuses the exp@rq check's own poly)

# ==================== @rq extension: log, sin/cos/tan, atan, asin/acos, composed ====================
#  Two different provenances, mirroring the @rh section's split:
#   - log's atanh-series coeffs and sin/cos's Taylor coeffs are EXACT closed
#     forms (1/(2k+3), (-1)^k/k!) -- confirmed to reproduce math.hoon's baked
#     hex bit-for-bit by regenerating them here (see the module's own
#     spot-check during development); using the formula is simpler than
#     transcribing 16-23 128-bit hex constants by hand.
#   - atan's and asin/acos's polynomials are genuine minimax FITS (fdlibm has
#     no quad-precision table to copy, unlike @rd/@rs), so an independent
#     re-fit is not guaranteed to reproduce the exact shipped bits.  These are
#     instead EXTRACTED PROGRAMMATICALLY from math.hoon's literal hex arrays
#     (see _extract_hoon_hex_list below) -- avoids hand-transcribing ~60
#     128-bit constants, which is exactly the kind of single-hex-nibble typo
#     that would silently corrupt a "measured ULP" number in this report.
#  Grid-based (not bit-pattern-exhaustive, matching check_exp_rq's own
#  convention): @rq has no practical exhaustive test (2^128 patterns), and
#  this file doesn't model @rq's raw NaN/Inf bit encoding at all (unlike the
#  @rh section above) -- consistent with how the existing check_exp_rq()
#  already tests a plain numeric grid, not decoded bit patterns.
import re as _re, os as _os
def _extract_hoon_hex_list(door_marker, arm_marker, end_marker, count):
    path = _os.path.join(_os.path.dirname(__file__), '..', 'desk', 'lib', 'math.hoon')
    text = open(path).read()
    i = text.index(door_marker); j = text.index(arm_marker, i); k = text.index(end_marker, j)
    toks = _re.findall(r'0x[0-9a-fA-F.]+', text[j:k])
    vals = [int(t[2:].replace('.', ''), 16) for t in toks]
    assert len(vals) == count, f"{door_marker}/{arm_marker}: expected {count} hex constants, found {len(vals)}"
    return vals

def qneg(x): return -mp.mpf(x)
def qfma(a, b, c): return qr(mp.mpf(a) * mp.mpf(b) + mp.mpf(c))   # extra working precision (dps=60 >> 113 bits) makes this safely single-rounded, matching a real fma
def q_head(x, keep_bits):
    """Zero out the low (112-keep_bits) mantissa bits of a 113-bit value
    (mirrors Hoon's `(dis s <mask>)` head/tail split used by +asin/+acos)."""
    if x == 0: return mp.mpf(0)
    m, e = mp.frexp(x)
    sig = int(mp.nint(m * mp.mpf(2) ** 113))
    masked = (sig >> (113 - keep_bits)) << (113 - keep_bits)
    return qr(mp.mpf(masked) * mp.mpf(2) ** (e - 113))

# ---- log/lr @rq (x=2^e*m + degree-22 EXACT atanh series; math.hoon ++log/++lr:rq) ----
LOG_COEFFS_Q = [qr(mp.mpf(1) / (2 * k + 3)) for k in range(23)]     # confirmed bit-identical to math.hoon's baked hex
SQRT2_Q = qr(mp.sqrt(2))
def lr_q(x):                           # x: qr()'d mpf, finite x>0 -> (e, l1)
    m, e = mp.frexp(x)                 # x = m*2^e, m in [0.5,1)
    m = qr(m * 2); e -= 1              # renormalize: m in [1,2)  (exact: doubling loses no bits)
    if m >= SQRT2_Q: m = qmul(m, mp.mpf('0.5')); e += 1
    f = qsub(m, mp.mpf(1))
    s = qdiv(f, qadd(m, mp.mpf(1)))
    z = qmul(s, s)
    p2 = qhorner(LOG_COEFFS_Q, z)
    r = qmul(qadd(z, z), p2)
    l1 = qsub(f, qmul(s, qsub(f, r)))
    return mp.mpf(e), l1
def log_q(x):
    x = qr(x)
    ef, l1 = lr_q(x)
    hi = qmul(ef, QLN2HI); lo = qmul(ef, QLN2LO)
    return qadd(hi, qadd(l1, lo))
LOG10_2_Q = qofhex(0x3ffd34413509f79fef311f12b35816f9)
INVLN10_Q = qofhex(0x3ffdbcb7b1526e50e32a6ab7555f5a68)
def log2_q(x):
    ef, lm = lr_q(qr(x))
    return qadd(ef, qmul(lm, QLOG2E))
def log10_q(x):
    ef, lm = lr_q(qr(x))
    return qadd(qmul(ef, LOG10_2_Q), qmul(lm, INVLN10_Q))

def check_log_rq():
    print("# log @rq: x=2^e*m reduction + degree-22 EXACT atanh series (formula-derived; confirmed vs math.hoon hex)")
    worst = mp.mpf(0); xw = None; n = 0
    for t in range(1, 200001, 37):
        x = qr(mp.mpf(t) / 100); g = log_q(x); tr = mp.log(x)
        ulp = mp.mpf(2) ** (mp.frexp(g)[1] - QP) if g != 0 else mp.mpf(2) ** -QP
        e = abs((g - tr) / ulp); n += 1
        if e > worst: worst, xw = e, x
    print(f"  max {float(worst):.3f} ULP at x={xw}  (sampled, n={n}, x in (0,2000])")
def check_log2log10_rq():
    print("# log-2/log-10 @rq: +lr split + linear combine (composed)")
    for name, fn, truefn in [('log-2', log2_q, lambda xm: mp.log(xm) / mp.log(2)),
                              ('log-10', log10_q, lambda xm: mp.log(xm) / mp.log(10))]:
        worst = mp.mpf(0); xw = None; n = 0
        for t in range(1, 200001, 37):
            x = qr(mp.mpf(t) / 100); g = fn(x); tr = truefn(x)
            ulp = mp.mpf(2) ** (mp.frexp(g)[1] - QP) if g != 0 else mp.mpf(2) ** -QP
            e = abs((g - tr) / ulp); n += 1
            if e > worst: worst, xw = e, x
        print(f"  {name}: max {float(worst):.3f} ULP at x={xw}  (sampled, n={n})")

# ---- sin/cos/tan @rq (quarter-turn reduction + degree-16 EXACT Taylor kernels; math.hoon ++rq-trig) ----
SC_Q = [qr(mp.mpf((-1) ** (k + 1)) / mp.factorial(2 * k + 3)) for k in range(16)]   # confirmed vs math.hoon hex
CC_Q = [qr(mp.mpf((-1) ** k) / mp.factorial(2 * k + 4)) for k in range(16)]         # confirmed vs math.hoon hex
INVPIO2_Q = qofhex(0x3ffe45f306dc9c882a53f84eafa3ea6a)
PIO2_1_Q  = qofhex(0x3fff921fb54442d18460000000000000)
PIO2_1T_Q = qofhex(0x3fc2313198a2e03707344a409382229a)
def ksin_q(xx, yy):
    z = qmul(xx, xx)
    r = qhorner(SC_Q[1:], z)
    v = qmul(z, xx)
    aa = qsub(qmul(mp.mpf('0.5'), yy), qmul(v, r))
    bb = qsub(qmul(z, aa), yy)
    dd = qsub(bb, qmul(v, SC_Q[0]))
    return qsub(xx, dd)
def kcos_q(xx, yy):
    z = qmul(xx, xx)
    rc = qhorner(CC_Q, z)
    hz = qmul(mp.mpf('0.5'), z)
    w2 = qsub(mp.mpf(1), hz)
    aa = qsub(qsub(mp.mpf(1), w2), hz)
    bb = qsub(qmul(qmul(z, z), rc), qmul(xx, yy))
    return qadd(w2, qadd(aa, bb))
def trig_fin_q(is_sin, ax, sb):
    q = int(mp.floor(qr(qmul(ax, INVPIO2_Q)) + mp.mpf('0.5')))
    qf = mp.mpf(abs(q))
    t = qsub(ax, qmul(qf, PIO2_1_Q))
    w = qmul(qf, PIO2_1T_Q)
    rhi = qsub(t, w); rlo = qsub(qsub(t, rhi), w)
    m = abs(q) & 3
    ks = ksin_q(rhi, rlo); kc = kcos_q(rhi, rlo)
    if is_sin:
        v = ks if m == 0 else kc if m == 1 else (qneg(ks) if m == 2 else qneg(kc))
        return qneg(v) if sb == 1 else v
    return kc if m == 0 else (qneg(ks) if m == 1 else (qneg(kc) if m == 2 else ks))
def sin_q(x):
    x = qr(x)
    if x == 0: return mp.mpf(0)
    neg = x < 0; ax = qneg(x) if neg else x
    return trig_fin_q(True, ax, 1 if neg else 0)
def cos_q(x):
    x = qr(x)
    ax = qneg(x) if x < 0 else x
    return trig_fin_q(False, ax, 0)
def tan_q(x): return qdiv(sin_q(x), cos_q(x))

def check_trig_rq(which):
    fn = sin_q if which == 'sin' else cos_q
    tru = mp.sin if which == 'sin' else mp.cos
    print(f"# {which} @rq: quarter-turn reduction + degree-16 EXACT Taylor kernel (formula-derived; confirmed vs math.hoon hex)")
    worst = mp.mpf(0); xw = None; n = 0
    for t in range(-200000, 200001, 37):
        x = qr(mp.mpf(t) / 100); g = fn(x); tr = tru(x)
        if abs(tr) < mp.mpf('1e-6'): continue
        ulp = mp.mpf(2) ** (mp.frexp(g)[1] - QP) if g != 0 else mp.mpf(2) ** -QP
        e = abs((g - tr) / ulp); n += 1
        if e > worst: worst, xw = e, x
    print(f"  max {float(worst):.3f} ULP at x={xw}  (sampled, n={n}, x in [-2000,2000])")
def check_tan_rq():
    print("# tan @rq: sin/cos ratio (formula-derived kernels)")
    worst = mp.mpf(0); xw = None; n = 0
    for t in range(-200000, 200001, 37):
        x = qr(mp.mpf(t) / 100); tr = mp.tan(x)
        if abs(tr) < mp.mpf('1e-6') or abs(tr) > mp.mpf('1e15'): continue
        g = tan_q(x)
        ulp = mp.mpf(2) ** (mp.frexp(g)[1] - QP) if g != 0 else mp.mpf(2) ** -QP
        e = abs((g - tr) / ulp); n += 1
        if e > worst: worst, xw = e, x
    print(f"  max {float(worst):.3f} ULP at x={xw}  (sampled, n={n})")

# ---- atan @rq (fdlibm breakpoint reduction + degree-30 minimax; math.hoon ++atan/++rq-atan) ----
AT_Q = [qofhex(h) for h in _extract_hoon_hex_list('++  rq-atan', '++  at\n', '++  atred', 31)]
ATB_Q = [qofhex(h) for h in (0x3ffdc000000000000000000000000000, 0x3ffe6000000000000000000000000000,
                              0x3fff3000000000000000000000000000, 0x40003800000000000000000000000000)]
ATHI_Q = [qofhex(h) for h in (0x3ffddac670561bb4f68adfc88bd97875,
                               0x3ffe921fb54442d18469898cc51701b8,
                               0x3ffef730bd281f69b200f10f5e197794,
                               0x3fff921fb54442d18469898cc51701b8)]
ATLO_Q = [qofhex(h) for h in (0x3f89a06dc282b0e4c39be01c59e2dcdd,
                               0x3f8bcd129024e088a67cc74020bbea64,
                               0xbf8bebe566c99ada9f231bccae27916c,
                               0x3f8ccd129024e088a67cc74020bbea64)]
def atan_ker_q(ax):
    if ax < ATB_Q[0]:
        xr, hi, lo, direct = ax, mp.mpf(0), mp.mpf(0), True
    elif ax < ATB_Q[1]:
        xr = qdiv(qsub(qadd(ax, ax), mp.mpf(1)), qadd(mp.mpf(2), ax)); hi, lo, direct = ATHI_Q[0], ATLO_Q[0], False
    elif ax < ATB_Q[2]:
        xr = qdiv(qsub(ax, mp.mpf(1)), qadd(ax, mp.mpf(1))); hi, lo, direct = ATHI_Q[1], ATLO_Q[1], False
    elif ax < ATB_Q[3]:
        xr = qdiv(qsub(ax, mp.mpf('1.5')), qadd(mp.mpf(1), qmul(mp.mpf('1.5'), ax))); hi, lo, direct = ATHI_Q[2], ATLO_Q[2], False
    else:
        xr = qneg(qdiv(mp.mpf(1), ax)); hi, lo, direct = ATHI_Q[3], ATLO_Q[3], False
    z = qmul(xr, xr)
    s = qmul(z, qhorner(AT_Q, z))
    if direct: return qsub(xr, qmul(xr, s))
    return qsub(hi, qsub(qsub(qmul(xr, s), lo), xr))
def atan_q(x):
    x = qr(x)
    if x == 0: return mp.mpf(0)
    neg = x < 0; ax = qneg(x) if neg else x
    r = atan_ker_q(ax)
    return qneg(r) if neg else r

def check_atan_rq():
    print("# atan @rq: fdlibm breakpoint reduction + degree-30 minimax (hex extracted from math.hoon)")
    worst = mp.mpf(0); xw = None; n = 0
    for t in range(-500000, 500001, 37):
        x = qr(mp.mpf(t) / 100); g = atan_q(x); tr = mp.atan(x)
        if abs(tr) < mp.mpf('1e-8'): continue
        ulp = mp.mpf(2) ** (mp.frexp(g)[1] - QP) if g != 0 else mp.mpf(2) ** -QP
        e = abs((g - tr) / ulp); n += 1
        if e > worst: worst, xw = e, x
    print(f"  max {float(worst):.3f} ULP at x={xw}  (sampled, n={n}, x in [-5000,5000])")

# ---- asin/acos @rq (fdlibm-style degree-30 rational-fit kernel; math.hoon ++asin/++acos/++rq-ainv) ----
RR_Q = [qofhex(h) for h in _extract_hoon_hex_list('++  rq-ainv', '++  rr\n', '++  asn', 31)]
PIO2_Q = qofhex(0x3fff921fb54442d18469898cc51701b8)
PIO2_LO_Q = qofhex(0x3f8ccd129024e088a67cc74020bbea64)
PIO4_Q = qofhex(0x3ffe921fb54442d18469898cc51701b8)
PI_Q = qofhex(0x4000921fb54442d18469898cc51701b8)
ASIN_THRESH_Q = qofhex(0x3ffef333333333333333333333333333)   # ~0.975
def rr_q(t): return qhorner(RR_Q, t)
def asin_q(x):
    x = qr(x)
    sgn = x < 0; ax = qneg(x) if sgn else x
    one = mp.mpf(1)
    if ax > one: return mp.mpf('nan')
    if ax == one:
        return qadd(qmul(x, PIO2_Q), qmul(x, PIO2_LO_Q))
    if ax < mp.mpf('0.5'):
        t = qmul(x, x)
        return qadd(x, qmul(x, rr_q(t)))
    w = qsub(one, ax); t = qmul(w, mp.mpf('0.5'))
    r = rr_q(t); s = qr(mp.sqrt(t))
    if ax >= ASIN_THRESH_Q:
        res = qsub(PIO2_Q, qsub(qmul(mp.mpf(2), qadd(s, qmul(s, r))), PIO2_LO_Q))
        return qneg(res) if sgn else res
    df = q_head(s, 56)
    c = qdiv(qsub(t, qmul(df, df)), qadd(s, df))
    p2 = qsub(qmul(mp.mpf(2), qmul(s, r)), qsub(PIO2_LO_Q, qmul(mp.mpf(2), c)))
    q2 = qsub(PIO4_Q, qmul(mp.mpf(2), df))
    res = qsub(PIO4_Q, qsub(p2, q2))
    return qneg(res) if sgn else res
def acos_q(x):
    x = qr(x)
    neg = x < 0; ax = qneg(x) if neg else x
    one = mp.mpf(1)
    if ax > one: return mp.mpf('nan')
    if ax == one:
        return mp.mpf(0) if not neg else qadd(PI_Q, qmul(mp.mpf(2), PIO2_LO_Q))
    if ax < mp.mpf('0.5'):
        z = qmul(x, x); r = rr_q(z)
        return qsub(PIO2_Q, qsub(x, qsub(PIO2_LO_Q, qmul(x, r))))
    if neg:
        z = qmul(qadd(one, x), mp.mpf('0.5')); s = qr(mp.sqrt(z)); r = rr_q(z)
        w = qsub(qmul(r, s), PIO2_LO_Q)
        return qsub(PI_Q, qmul(mp.mpf(2), qadd(s, w)))
    z = qmul(qsub(one, x), mp.mpf('0.5')); s = qr(mp.sqrt(z))
    df = q_head(s, 56)
    c = qdiv(qsub(z, qmul(df, df)), qadd(s, df))
    r = rr_q(z); w = qadd(qmul(r, s), c)
    return qmul(mp.mpf(2), qadd(df, w))

def check_ainv_rq(which):
    fn = asin_q if which == 'asin' else acos_q
    tru = mp.asin if which == 'asin' else mp.acos
    print(f"# {which} @rq: fdlibm-style degree-30 rational-fit kernel (hex extracted from math.hoon)")
    worst = mp.mpf(0); xw = None; n = 0
    for t in range(-1000000, 1000001, 71):
        x = qr(mp.mpf(t) / 1000000); g = fn(x); tr = tru(x)
        if abs(tr) < mp.mpf('1e-9'): continue
        ulp = mp.mpf(2) ** (mp.frexp(g)[1] - QP) if g != 0 else mp.mpf(2) ** -QP
        e = abs((g - tr) / ulp); n += 1
        if e > worst: worst, xw = e, x
    print(f"  max {float(worst):.3f} ULP at x={xw}  (sampled, n={n}, x in [-1,1])")

# ---- composed @rq: sqt/cbrt/pow/pow-n/atan2 ----
def sqt_q(x):
    x = qr(x)
    if x == 0: return x
    g = qr(mp.sqrt(x))
    h = qdiv(mp.mpf('0.5'), g)
    r = qfma(qneg(g), g, x)
    return qfma(h, r, g)
THIRD_Q = qofhex(0x3ffd5555555555555555555555555555)
def cbt_q(x):
    x = qr(x)
    if x == 0: return x
    ax = -x if x < 0 else x
    r = exp_q(qmul(log_q(ax), THIRD_Q), EXP_P_Q)
    return qneg(r) if x < 0 else r
def pow_n_q(x, n):
    n = int(n)
    if n == 0: return mp.mpf(1)
    p = qr(x)
    for _ in range(n - 1): p = qmul(p, x)
    return p
def pow_q(x, n):
    if n == qr(mp.nint(n)) and n > 0:
        return pow_n_q(x, int(n))
    return exp_q(qmul(n, log_q(x)), EXP_P_Q)
def atan2_q(y, x):
    y, x = qr(y), qr(x)
    if x > 0: return atan_q(qdiv(y, x))
    if x < 0 and y >= 0: return qadd(atan_q(qdiv(y, x)), PI_Q)
    if x < 0 and y < 0: return qsub(atan_q(qdiv(y, x)), PI_Q)
    if x == 0 and y > 0: return PIO2_Q
    if x == 0 and y < 0: return qneg(PIO2_Q)
    return mp.mpf(0)

def _qulp_sweep(fn, truefn, xs, skip_small=mp.mpf(0)):
    worst = mp.mpf(0); xw = None; n = 0
    for x in xs:
        g = fn(x)
        if g != g or not mp.isfinite(g): continue
        tr = truefn(x)
        if not mp.isfinite(tr) or abs(tr) < skip_small: continue
        ulp = mp.mpf(2) ** (mp.frexp(g)[1] - QP) if g != 0 else mp.mpf(2) ** -QP
        e = abs((g - tr) / ulp); n += 1
        if e > worst: worst, xw = e, x
    return worst, xw, n

def check_sqt_rq():
    print("# sqt @rq: Markstein-corrected sqrt (composed)")
    xs = [qr(mp.mpf(t) / 1000) for t in range(1, 2000001, 977)]
    worst, xw, n = _qulp_sweep(sqt_q, mp.sqrt, xs)
    print(f"  max {float(worst):.3f} ULP at x={xw}  (sampled, n={n}, x in (0,2000])")
def check_cbrt_rq():
    print("# cbrt @rq: sign(x)*exp(log|x|/3) (composed)")
    xs = [qr(mp.mpf(t) / 1000) for t in range(-2000000, 2000001, 977) if t != 0]
    def _cbrt_true_q(xm): return -mp.cbrt(-xm) if xm < 0 else mp.cbrt(xm)
    worst, xw, n = _qulp_sweep(cbt_q, _cbrt_true_q, xs)
    print(f"  max {float(worst):.3f} ULP at x={xw}  (sampled, n={n}, x in [-2000,2000])")
def check_pow_rq():
    print("# pow/pow-n @rq: pow-n repeated-multiply (int n>0) else exp(n*log x) (composed)")
    worst = mp.mpf(0); xw = None; n = 0
    for xt in range(1, 401, 3):
        x = qr(mp.mpf(xt) / 20)
        for nt in range(-400, 401, 11):
            nn = qr(mp.mpf(nt) / 20)
            g = pow_q(x, nn)
            if g != g or not mp.isfinite(g): continue
            tr = x ** nn
            if not mp.isfinite(tr) or tr == 0: continue
            ulp = mp.mpf(2) ** (mp.frexp(g)[1] - QP) if g != 0 else mp.mpf(2) ** -QP
            e = abs((g - tr) / ulp); n += 1
            if e > worst: worst, xw = e, (x, nn)
    print(f"  max {float(worst):.3f} ULP at (x,n)={xw}  (sampled grid, n={n}, x in (0,20], n in [-20,20])")
def check_atan2_rq():
    print("# atan2 @rq: quadrant dispatch over +atan (composed)")
    worst = mp.mpf(0); xw = None; n = 0
    for yt in range(-200, 201, 7):
        for xt in range(-200, 201, 7):
            y = qr(mp.mpf(yt) / 10); x = qr(mp.mpf(xt) / 10)
            if y == 0 and x == 0: continue
            g = atan2_q(y, x)
            tr = mp.atan2(y, x)
            ulp = mp.mpf(2) ** (mp.frexp(g)[1] - QP) if g != 0 else mp.mpf(2) ** -QP
            e = abs((g - tr) / ulp); n += 1
            if e > worst: worst, xw = e, (y, x)
    print(f"  max {float(worst):.3f} ULP at (y,x)={xw}  (sampled grid, n={n}, [-20,20]x[-20,20])")

if __name__ == '__main__':
    fn = sys.argv[1] if len(sys.argv) > 1 else 'exp'
    {'exp': check_exp, 'exp-rs': check_exp_rs,
     'log': check_log, 'log-rs': check_log_rs,
     'sin': lambda: check_trig('sin'), 'cos': lambda: check_trig('cos'),
     'sin-rs': lambda: check_trig_rs('sin'), 'cos-rs': lambda: check_trig_rs('cos'),
     'atan': check_atan, 'atan-rs': check_atan_rs, 'asin': check_asin, 'acos': check_acos,
     'asin-rs': lambda: check_ainv_rs('asin'), 'acos-rs': lambda: check_ainv_rs('acos'),
     'tan': check_tan, 'tan-rs': check_tan_rs, 'exp-rq': check_exp_rq,
     # -- new: @rd/@rs composed ops --
     'log2log10-rd': check_log2log10_rd, 'log2log10-rs': check_log2log10_rs,
     'sqt-rd': check_sqt_rd, 'sqt-rs': check_sqt_rs,
     'cbrt-rd': check_cbrt_rd, 'cbrt-rs': check_cbrt_rs,
     'pow-rd': check_pow_rd, 'pow-rs': check_pow_rs,
     'atan2-rd': check_atan2_rd, 'atan2-rs': check_atan2_rs,
     # -- new: @rh (exhaustive over 2^16) --
     'exp-rh': check_exp_rh, 'log-rh': check_log_rh, 'log2log10-rh': check_log2log10_rh,
     'sin-rh': lambda: check_trig_rh('sin'), 'cos-rh': lambda: check_trig_rh('cos'), 'tan-rh': check_tan_rh,
     'atan-rh': check_atan_rh, 'asin-rh': lambda: check_ainv_rh('asin'), 'acos-rh': lambda: check_ainv_rh('acos'),
     'sqt-rh': check_sqt_rh, 'cbrt-rh': check_cbrt_rh, 'pow-rh': check_pow_rh, 'atan2-rh': check_atan2_rh,
     # -- new: @rq extension --
     'log-rq': check_log_rq, 'log2log10-rq': check_log2log10_rq,
     'sin-rq': lambda: check_trig_rq('sin'), 'cos-rq': lambda: check_trig_rq('cos'), 'tan-rq': check_tan_rq,
     'atan-rq': check_atan_rq, 'asin-rq': lambda: check_ainv_rq('asin'), 'acos-rq': lambda: check_ainv_rq('acos'),
     'sqt-rq': check_sqt_rq, 'cbrt-rq': check_cbrt_rq, 'pow-rq': check_pow_rq, 'atan2-rq': check_atan2_rq,
     }[fn]()
