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

if __name__ == '__main__':
    fn = sys.argv[1] if len(sys.argv) > 1 else 'exp'
    {'exp': check_exp, 'exp-rs': check_exp_rs,
     'log': check_log, 'log-rs': check_log_rs}[fn]()
