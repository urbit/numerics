#!/usr/bin/env python3
"""Offline oracle harness for /lib/fixed (fixed-point Q a.b).

Cross-checks a Python mirror of the Hoon /lib/fixed algorithm against TWO
independent oracles:

  1. exact rationals (`fractions.Fraction`) -- a fixed-point value is exactly
     scaledval / 2^b, so this is ground truth with no rounding of its own.
  2. the `spfpm` library (the `FixedPoint` module), an established third-party
     fixed-point implementation.

A fixed-point number here is a two's-complement integer of N = a+b+1 bits
interpreted as value/2^b (sign bit at N-1).  Arithmetic is signed and modular
(wraps mod 2^N), matching the Hoon library and /lib/twoc.

    pip install spfpm
    python3 libmath/tools/fixed_check.py
"""
from fractions import Fraction
try:
    from FixedPoint import FXfamily, FXnum
except ImportError:
    FXnum = None

# ---- mirror of the Hoon /lib/fixed (delegates to twoc width N=a+b+1) ----

def nof(a, b):              return a + b + 1
def s2t(s, n):             return s % (1 << n)              # signed int -> N bits
def t2s(t, n):             return t - (1 << n) if t >> (n - 1) else t

def f_add(x, xp, y, yp):
    n = nof(*xp); return (x + y) % (1 << n)
def f_sub(x, xp, y, yp):
    n = nof(*xp); return (x - y) % (1 << n)
def f_neg(x, xp):
    n = nof(*xp); return (-t2s(x, n)) % (1 << n)
def f_mul(x, xp, y, yp):
    pa, pb = xp[0] + yp[0] + 1, xp[1] + yp[1]; n = nof(pa, pb)
    return ((t2s(x, nof(*xp)) * t2s(y, nof(*yp))) % (1 << n), (pa, pb))
def f_div(x, xp, y, yp):
    n = nof(*xp); sx = t2s(x, n); sy = t2s(y, nof(*yp))
    num = sx << xp[1]
    q = abs(num) // abs(sy)
    if (num < 0) != (sy < 0): q = -q
    return (q % (1 << n), xp)
def f_scale(x, xp, yp):
    nx, ny = nof(*xp), nof(*yp); sx = t2s(x, nx)
    sg = sx >= 0; mg = abs(sx)
    mg2 = mg << (yp[1] - xp[1]) if yp[1] > xp[1] else mg >> (xp[1] - yp[1])
    v = mg2 if sg else -mg2
    return v % (1 << ny)

# ---- value views ----

def fxval(bits, a, b):     return Fraction(t2s(bits, nof(a, b)), 1 << b)  # exact value

# ---- checks ----

def near(x, y, tol=Fraction(0)):
    return abs(x - y) <= tol

def check_exact():
    """Mirror ops vs exact-rational expectation.

    The reference accounts for the library's MODULAR wrap: a result is the
    true rational re-encoded into the result width and decoded back, so an
    out-of-range value wraps exactly as the Hoon does.
    """
    a, b = 8, 8
    samples = [Fraction(3, 2), Fraction(-3, 2), Fraction(9, 4), Fraction(-1),
               Fraction(5), Fraction(-7, 1), Fraction(1, 256)]
    def enc(v, bb=b):  return s2t(int(v * (1 << bb)), nof(a, b))   # value -> q8.8 bits
    def wrap(v, ra, rb):                       # true value -> (re-encoded) value at qRa.Rb
        return fxval(s2t(int(v * (1 << rb)), nof(ra, rb)), ra, rb)
    def trunc(q):                              # truncate toward zero to 1/2^b
        return Fraction(int(q * 256), 256) if q >= 0 else -Fraction(int(-q * 256), 256)
    ok = True
    for x in samples:
        for y in samples:
            bx, by = enc(x), enc(y)
            if fxval(f_add(bx, (a, b), by, (a, b)), a, b) != wrap(x + y, a, b):
                ok = False; print("add", x, y)
            if fxval(f_sub(bx, (a, b), by, (a, b)), a, b) != wrap(x - y, a, b):
                ok = False; print("sub", x, y)
            mb, mp = f_mul(bx, (a, b), by, (a, b))
            if fxval(mb, mp[0], mp[1]) != wrap(x * y, mp[0], mp[1]):
                ok = False; print("mul", x, y)
            if y != 0:
                db, dp = f_div(bx, (a, b), by, (a, b))
                if fxval(db, dp[0], dp[1]) != wrap(trunc(x / y), dp[0], dp[1]):
                    ok = False; print("div", x, y, float(fxval(db, dp[0], dp[1])))
    print("exact-rational: add/sub/mul/div over samples^2 (with wrap):",
          "OK" if ok else "MISMATCH")
    return ok

def check_spfpm():
    """Parity with the spfpm library on IN-RANGE results.  spfpm is a CHECKED
    (overflow-raising) type, so it can only witness cases that don't wrap;
    we skip any op whose result leaves q8.8 range and let check_exact cover
    the wrapping cases.
    """
    if FXnum is None:
        print("spfpm: SKIP (not installed)"); return True
    a, b = 8, 8
    fam = FXfamily(b, a)
    # q8.8 holds signed scaledval in [-2^15, 2^15-1].  spfpm's FXfamily(8,8)
    # validates the same magnitude range, so we keep operands and results
    # strictly inside it (it raises on overflow rather than wrapping).
    lo, hi = -(1 << 15), (1 << 15) - 1
    import random; random.seed(3); ok = True; checked = 0
    def mk(s):
        v = FXnum(0, fam); v.scaledval = s; return v
    while checked < 8000:
        sx = random.randrange(lo, hi + 1)
        sy = random.randrange(lo, hi + 1)
        bx, by = s2t(sx, nof(a, b)), s2t(sy, nof(a, b))
        for op, res in (('add', sx + sy), ('sub', sx - sy)):
            if not (lo <= res <= hi):  continue
            r = (mk(sx) + mk(sy)) if op == 'add' else (mk(sx) - mk(sy))
            mine = f_add(bx, (a, b), by, (a, b)) if op == 'add' else f_sub(bx, (a, b), by, (a, b))
            if fxval(mine, a, b) != Fraction(r.scaledval, 1 << b):
                ok = False; print(op, sx, sy)
            checked += 1
    print(f"spfpm: add/sub parity on {checked} in-range q8.8 cases:",
          "OK" if ok else "MISMATCH")
    return ok

if __name__ == '__main__':
    e = check_exact()
    s = check_spfpm()
    print("ALL PASS:", e and s)
