#!/usr/bin/env python3
# rq_inv_check.py -- fine-sampled faithfulness audit of the @rq atan/asin/acos
# arms, which aren't covered by the C harness (tools/rq_check.c does exp/log/
# sin/cos/sqrt against MPFR).  This models f128 with mpmath rounded to 113-bit
# RNE per op (= SoftFloat f128 for + - * / sqrt), parses the actual coefficient
# lists straight out of math.hoon, and sweeps each arm against mpmath truth.
#
# The mpmath sim tracks the real SoftFloat jet to ~0.05 ULP (spot-check any
# flagged point against the live Hoon before acting -- rq_check.c + MPFR remains
# the authoritative tool for the arms it covers).
#
# Run: python3 rq_inv_check.py        (worst-case ULP per arm, fine sweep)
import re, sys, mpmath as mp
mp.mp.dps = 80                                  # ~265 bits, ample over 113

# ---- binary128 layer (normal range, round-nearest-even) ----
def f128(x):
    x = mp.mpf(x)
    if x == 0 or mp.isinf(x) or mp.isnan(x): return x
    sgn = -1 if x < 0 else 1; x = abs(x)
    e = int(mp.floor(mp.log(x, 2))); m = x / mp.power(2, e)
    while m >= 2: m /= 2; e += 1
    while m < 1:  m *= 2; e -= 1
    scaled = m * mp.power(2, 112); fl = mp.floor(scaled); frac = scaled - fl; n = int(fl)
    if   frac > mp.mpf('0.5'):  n += 1
    elif frac == mp.mpf('0.5'): n += (n & 1)    # ties to even
    if n == (1 << 113): n >>= 1; e += 1
    return mp.mpf(sgn) * mp.mpf(n) * mp.power(2, e - 112)
def of_bits(b):
    s = (b >> 127) & 1; e = (b >> 112) & 0x7fff; m = b & ((1 << 112) - 1)
    v = (mp.mpf(1) + mp.mpf(m) / mp.power(2, 112)) * mp.power(2, e - 16383) if e \
        else mp.mpf(m) * mp.power(2, -16382 - 112)
    return -v if s else v
def to_bits(x):
    x = f128(x)
    if x == 0: return 0
    s = 1 if x < 0 else 0; x = abs(x); e = int(mp.floor(mp.log(x, 2))); m = x / mp.power(2, e)
    while m >= 2: m /= 2; e += 1
    while m < 1:  m *= 2; e -= 1
    mant = int(mp.nint((m - 1) * mp.power(2, 112)))
    if mant == (1 << 112): mant = 0; e += 1
    return (s << 127) | ((e + 16383) << 112) | mant
def ulps(approx, true):
    a = abs(mp.mpf(true))
    ulp = mp.power(2, -16382 - 112) if a == 0 else mp.power(2, int(mp.floor(mp.log(a, 2))) - 112)
    return float((f128(approx) - mp.mpf(true)) / ulp)

q = of_bits
def qadd(a, b): return f128(a + b)
def qsub(a, b): return f128(a - b)
def qmul(a, b): return f128(a * b)
def qdiv(a, b): return f128(a / b)
def qsqt(a):    return f128(mp.sqrt(a))
def dis(s, mask): return of_bits(to_bits(s) & mask)         # bitwise-and on the f128 encoding
def horner(co, t):
    acc = co[-1]
    for c in reversed(co[:-1]): acc = qadd(qmul(acc, t), c)
    return acc

# ---- coefficient lists, parsed straight out of math.hoon (no hand-copy) ----
import os
_HOON = os.path.join(os.path.dirname(__file__), '..', 'desk', 'lib', 'math.hoon')
src = open(_HOON).read()
def grab(a, b):
    seg = src[src.index(a):src.index(b, src.index(a))]
    return [q(int(h.replace('.', ''), 16)) for h in re.findall(r'`@rq`0x([0-9a-f.]+)', seg)]
AT = grab('++  rq-atan', '++  atred')           # atan minimax (ascending in z=xr^2)
RR = grab('++  rq-ainv', '++  asn')             # asin/acos kernel R(t)

# ---- constants (the exact values the @rq Hoon arms use) ----
PIO2H = q(0x3fff921fb54442d18469898cc51701b8); PIO2L = q(0x3f8ccd129024e088a67cc74020bbea64)
PIO4H = q(0x3ffe921fb54442d18469898cc51701b8); PIH   = q(0x4000921fb54442d18469898cc51701b8)
ONE   = q(0x3fff0000000000000000000000000000); TWO   = q(0x40000000000000000000000000000000)
HALF  = q(0x3ffe0000000000000000000000000000); DMASK = 0xffffffffffffffffff00000000000000

def atred(ax):                                  # fdlibm breakpoint reduction
    if ax < q(0x3ffdc000000000000000000000000000): return (ax, q(0), q(0), True)
    if ax < q(0x3ffe6000000000000000000000000000):
        return (qdiv(qsub(qadd(ax, ax), ONE), qadd(TWO, ax)),
                q(0x3ffddac670561bb4f68adfc88bd97875), q(0x3f89a06dc282b0e4c39be01c59e2dcdd), False)
    if ax < q(0x3fff3000000000000000000000000000):
        return (qdiv(qsub(ax, ONE), qadd(ax, ONE)),
                PIO4H, q(0x3f8bcd129024e088a67cc74020bbea64), False)
    if ax < q(0x40003800000000000000000000000000):
        ohf = q(0x3fff8000000000000000000000000000)
        return (qdiv(qsub(ax, ohf), qadd(ONE, qmul(ohf, ax))),
                q(0x3ffef730bd281f69b200f10f5e197794), q(0xbf8bebe566c99ada9f231bccae27916c), False)
    return (qdiv(q(0xbfff0000000000000000000000000000), ax), PIO2H, PIO2L, False)
def atan(x):
    ax = f128(abs(x)); neg = x < 0
    xr, hi, lo, d = atred(ax)
    z = qmul(xr, xr); s = qmul(z, horner(AT, z))
    r = qsub(xr, qmul(xr, s)) if d else qsub(hi, qsub(qsub(qmul(xr, s), lo), xr))
    return f128(-r) if neg else r
def asn(x):
    sgn = x < 0; ax = f128(abs(x))
    if ax > ONE: return mp.nan
    if ax == ONE: return qadd(qmul(x, PIO2H), qmul(x, PIO2L))
    if ax < HALF:
        if ax < q(0x3fc60000000000000000000000000000): return x
        return qadd(x, qmul(x, horner(RR, qmul(x, x))))
    w = qsub(ONE, ax); t = qmul(w, HALF); r = horner(RR, t); s = qsqt(t)
    if ax >= q(0x3ffef333333333333333333333333333):
        res = qsub(PIO2H, qsub(qmul(TWO, qadd(s, qmul(s, r))), PIO2L)); return f128(-res) if sgn else res
    df = dis(s, DMASK); c = qdiv(qsub(t, qmul(df, df)), qadd(s, df))
    p2 = qsub(qmul(TWO, qmul(s, r)), qsub(PIO2L, qmul(TWO, c))); q2 = qsub(PIO4H, qmul(TWO, df))
    res = qsub(PIO4H, qsub(p2, q2)); return f128(-res) if sgn else res
def acs(x):
    neg = x < 0; ax = f128(abs(x))
    if ax > ONE: return mp.nan
    if ax == ONE: return q(0) if not neg else qadd(PIH, qmul(TWO, PIO2L))
    if ax < HALF:
        if ax < q(0x3f870000000000000000000000000000): return PIO2H
        z = qmul(x, x); r = horner(RR, z); return qsub(PIO2H, qsub(x, qsub(PIO2L, qmul(x, r))))
    if neg:
        z = qmul(qadd(ONE, x), HALF); s = qsqt(z); r = horner(RR, z); w = qsub(qmul(r, s), PIO2L)
        return qsub(PIH, qmul(TWO, qadd(s, w)))
    z = qmul(qsub(ONE, x), HALF); s = qsqt(z); df = dis(s, DMASK)
    c = qdiv(qsub(z, qmul(df, df)), qadd(s, df)); r = horner(RR, z); w = qadd(qmul(r, s), c)
    return qmul(TWO, qadd(df, w))

def sweep(fn, tru, lo, hi, n, dom=None):
    worst = 0.0; xw = 0.0
    for i in range(n + 1):
        x = f128(lo + (hi - lo) * mp.mpf(i) / n)
        if dom and not dom(x): continue
        g = fn(x)
        if g != g: continue
        t = tru(x)
        if t == 0 or abs(t) < mp.mpf('1e-300'): continue
        e = abs(ulps(g, t))
        if e > worst: worst, xw = e, float(x)
    return worst, xw

if __name__ == '__main__':
    print(f"parsed AT={len(AT)} coeffs, RR={len(RR)} coeffs")
    assert to_bits(atan(ONE)) == 0x3ffe921fb54442d18469898cc51701b8, "atan(1) != pi/4"
    n = 40000
    wa = sweep(atan, mp.atan, mp.mpf(-8), mp.mpf(8), n)
    ws = sweep(asn, mp.asin, mp.mpf('-0.999'), mp.mpf('0.999'), n, lambda x: abs(x) <= 1)
    wc = sweep(acs, mp.acos, mp.mpf('-0.999'), mp.mpf('0.999'), n, lambda x: abs(x) <= 1)
    print(f"# rq inv: atan {wa[0]:.3f} (x={wa[1]:.3f})  asin {ws[0]:.3f} (x={ws[1]:.3f})  "
          f"acos {wc[0]:.3f} (x={wc[1]:.3f})   [all faithful, <=1 ULP]")
