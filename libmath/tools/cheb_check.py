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

if __name__ == '__main__':
    fn = sys.argv[1] if len(sys.argv) > 1 else 'exp'
    {'exp': check_exp, 'exp-rs': check_exp_rs,
     'log': check_log, 'log-rs': check_log_rs,
     'sin': lambda: check_trig('sin'), 'cos': lambda: check_trig('cos'),
     'sin-rs': lambda: check_trig_rs('sin'), 'cos-rs': lambda: check_trig_rs('cos'),
     'atan': check_atan, 'atan-rs': check_atan_rs, 'asin': check_asin, 'acos': check_acos,
     'asin-rs': lambda: check_ainv_rs('asin'), 'acos-rs': lambda: check_ainv_rs('acos')}[fn]()
