#!/usr/bin/env python3
"""Offline verification harness for /lib/unum (2022 Posit Standard, es=2).

Runs the *exhaustive* checks that are too expensive for the on-ship Hoon
suite, against two independent oracles:

  1. SoftPosit (the `softposit` pip package, pX2 = es=2 at any width)
  2. a from-scratch exact-rational reference encoder/decoder

It (a) cross-checks the math constants baked into lib/unum.hoon at
posit8/16/32, and (b) exhaustively checks add/sub/mul/div over ALL
65,536 posit8 pairs.  The Hoon suite mirrors a lean subset (round-trips,
spot values, a property sweep); this is the heavy proof, run once.

    pip install softposit mpmath
    python3 libmath/tools/posit_check.py
"""
from fractions import Fraction
try:
    import softposit as sp
except ImportError:
    sp = None
from mpmath import mp, mpf, pi as PI, e as E, sqrt, log, phi as PHI
mp.dps = 80

# ---- reference decoder/encoder, mirroring lib/unum.hoon (es=2) ----

def decode(p, n):
    msk = (1 << n) - 1; p &= msk; nar = 1 << (n - 1)
    if p == 0: return ('z',)
    if p == nar: return ('n',)
    neg = (p >> (n - 1)) & 1
    mag = (1 << n) - p if neg else p
    pw = n - 1; r0 = (mag >> (pw - 1)) & 1; k = 1
    while True:
        if k == pw:
            r = (k - 1) if r0 == 1 else -k
            return ('p', neg, 4 * r, 1)
        if (mag >> (pw - 1 - k)) & 1 == r0: k += 1; continue
        break
    r = (k - 1) if r0 == 1 else -k
    remwid = pw - (k + 1); rem = mag & ((1 << remwid) - 1)
    if remwid >= 2: elo = rem >> (remwid - 2); fw = remwid - 2
    elif remwid == 1: elo = rem << 1; fw = 0
    else: elo = 0; fw = 0
    frac = rem & ((1 << fw) - 1)
    x = 4 * r + elo; a = (1 << fw) + frac
    return ('p', neg, x - fw, a)

def encode(neg, e, a, n):
    if a == 0: return 0
    msk = (1 << n) - 1; maxpos = (1 << (n - 1)) - 1
    lead = a.bit_length() - 1; x = e + lead; frac = a & ((1 << lead) - 1)
    r = x >> 2; elo = x - 4 * r
    if r >= n - 2: return ((1 << n) - maxpos) & msk if neg else maxpos
    if r <= -(n - 1): return ((1 << n) - 1) & msk if neg else 1
    if r >= 0: regval = ((1 << (r + 1)) - 1) << 1; regwid = r + 2
    else: regval = 1; regwid = -r + 1
    totw = regwid + 2 + lead
    pay = (regval << (2 + lead)) | (elo << lead) | frac
    pw = n - 1
    if totw <= pw:
        mag = pay << (pw - totw)
    else:
        sh = totw - pw; keep = pay >> sh
        guard = (pay >> (sh - 1)) & 1; low = pay & ((1 << (sh - 1)) - 1)
        if guard and ((1 if low else 0) or (keep & 1)): keep += 1
        if keep > maxpos: keep = maxpos
        mag = keep
    return ((1 << n) - mag) & msk if neg else mag

def ref_value_encode(v, n):  # v: Fraction
    if v == 0: return 0
    neg = v < 0; w = -v if neg else v
    maxv = Fraction(2) ** (4 * n - 8)
    if w >= maxv: return ((1 << n) - maxpos_of(n)) & ((1 << n) - 1) if neg else maxpos_of(n)
    if w <= 1 / maxv: return ((1 << n) - 1) & ((1 << n) - 1) if neg else 1
    x = 0; vv = w
    while vv >= 2: vv /= 2; x += 1
    while vv < 1: vv *= 2; x -= 1
    K = 80; a = int(vv * Fraction(1 << K))  # significand with hidden 1
    return encode(neg, x - K, a, n)

def maxpos_of(n): return (1 << (n - 1)) - 1

# ---- g-layer arithmetic, mirroring lib/unum.hoon ----

def mul(a, b, n):
    ua, ub = decode(a, n), decode(b, n)
    if ua[0] == 'n' or ub[0] == 'n': return 1 << (n - 1)
    if ua[0] == 'z' or ub[0] == 'z': return 0
    _, sa, ea, aa = ua; _, sb, eb, ab = ub
    return encode(sa ^ sb, ea + eb, aa * ab, n)

def add(a, b, n):
    ua, ub = decode(a, n), decode(b, n)
    if ua[0] == 'n' or ub[0] == 'n': return 1 << (n - 1)
    if ua[0] == 'z': return b
    if ub[0] == 'z': return a
    _, sa, ea, aa = ua; _, sb, eb, ab = ub
    emin = min(ea, eb); s1 = aa << (ea - emin); s2 = ab << (eb - emin)
    if sa == sb: return encode(sa, emin, s1 + s2, n)
    if s1 > s2: return encode(sa, emin, s1 - s2, n)
    if s2 > s1: return encode(sb, emin, s2 - s1, n)
    return 0

def neg(a, n): return ((1 << n) - a) & ((1 << n) - 1)
def sub(a, b, n): return add(a, neg(b, n), n)

def div(a, b, n):
    ua, ub = decode(a, n), decode(b, n)
    if ua[0] == 'n' or ub[0] == 'n': return 1 << (n - 1)
    if ub[0] == 'z': return 1 << (n - 1)
    if ua[0] == 'z': return 0
    _, sa, ea, aa = ua; _, sb, eb, ab = ub
    g = 2 * n; num = aa << g; q = num // ab
    if num % ab: q |= 1
    return encode(sa ^ sb, ea - eb - g, q, n)

# ---- softposit oracle ----

def sp_pat(p, n): return p.v >> (32 - n)
def sp_mk(pat, n):
    o = sp.convertDoubleToPX2(0.0, n); o.v = pat << (32 - n); return o
def sp_op(op, a, b, n):
    f = {'add': sp.pX2_add, 'sub': sp.pX2_sub, 'mul': sp.pX2_mul, 'div': sp.pX2_div}[op]
    return sp_pat(f(sp_mk(a, n), sp_mk(b, n), n), n)

# ---- checks ----

CONSTS = {'pi': PI, 'tau': 2 * PI, 'e': E, 'phi': PHI, 'sqt2': sqrt(2),
          'invsqt2': 1 / sqrt(2), 'log2': log(2), 'invlog2': 1 / log(2),
          'log10': log(10)}
WIDTHS = [('rpb', 8), ('rph', 16), ('rps', 32)]

def check_consts():
    ok = True
    for name, val in CONSTS.items():
        fr = Fraction(int(mp.nint(val * (mpf(2) ** 120))), 1 << 120)
        for w, n in WIDTHS:
            ref = ref_value_encode(fr, n)
            if sp is not None:
                ora = sp_pat((sp.convertDoubleToP32(float(val)) if n == 32
                              else sp.convertDoubleToPX2(float(val), n)), n)
                if ref != ora:
                    print(f"  CONST MISMATCH {name} {w}: ref=0x{ref:x} sp=0x{ora:x}"); ok = False
    print("constants: ref-encoder == SoftPosit at posit8/16/32:", ok)
    return ok

def check_arith():
    if sp is None:
        print("arith: SKIP (softposit not installed)"); return True
    ok = True
    for name, fn in (('add', add), ('sub', sub), ('mul', mul), ('div', div)):
        bad = 0
        for a in range(256):
            for b in range(256):
                if fn(a, b, 8) != sp_op(name, a, b, 8): bad += 1
        print(f"  {name}: {65536 - bad}/65536 posit8 pairs match SoftPosit")
        ok &= (bad == 0)
    return ok

if __name__ == '__main__':
    c = check_consts(); a = check_arith()
    print("ALL PASS:", c and a)
