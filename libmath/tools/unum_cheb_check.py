#!/usr/bin/env python3
"""Chebyshev/minimax transcendental harness for /lib/unum (posits), mirroring
`cheb_check.py`'s role for /lib/math.

Design (see numerics libmath/NEXT-STEPS.md, "Roadmap -- 2026-07-03
(Chebyshev transcendentals)"):
  - Unlike /lib/math's IEEE floats (fixed-width hardware registers, forcing a
    hi/lo Cody-Waite constant split), Hoon `@` atoms are arbitrary-precision.
    So this algorithm-of-record does exact `Fraction` arithmetic throughout
    range reduction and polynomial evaluation -- no intermediate rounding,
    just like a Hoon implementation that lets its fixed-point atoms grow and
    rounds ONLY ONCE at the end via `+bit`.  The only place finite precision
    enters is *storing* constants (ln2, pi/2, minimax coefficients) as fixed-
    point literals with WBITS bits of fraction -- generous enough (128 bits)
    that this quantization is far below posit32's worst-case ULP.
  - `decode`/`encode`/`ref_value_encode` are lib/unum.hoon's `+sea`/`+bit`,
    reused verbatim from `posit_check.py` (the existing faithful reference
    model) -- this file's new machinery is the range-reduction + minimax
    kernel that sits between them.
  - mpmath is the INDEPENDENT truth: it designs the coefficients (via
    mp.chebyfit on the reduced domain) and proves the faithful-rounding
    (<=1 ULP) bound.  It is NOT the test oracle; the algorithm-of-record
    (this file) is -- exactly cheb_check.py's division of labor.

Run:  python3 libmath/tools/unum_cheb_check.py exp
Requires: mpmath.
"""
import sys
from fractions import Fraction
import mpmath as mp
mp.mp.dps = 60

from posit_check import decode, encode, ref_value_encode, my_sqrt, maxpos_of
from posit_check import mul as pc_mul, sub as pc_sub, div as pc_div, neg as pc_neg

WIDTHS = [('rpb', 8), ('rph', 16), ('rps', 32)]
WBITS = 128                                    # fixed-point precision for baked constants

def quant_w(mv, w):
    """mpmath value -> Fraction, rounded to `w` bits of fraction (simulates a
    baked hex literal atom -- what the Hoon constant actually stores)."""
    return Fraction(int(mp.nint(mv * (mp.mpf(2) ** w))), 1 << w)

def raw_w(fr, w):
    """Fraction -> the integer numerator AS IF its denominator were exactly
    2^w (undoes Fraction's automatic gcd-reduction).  `Fraction` silently
    shrinks the denominator when the rounded numerator is even (e.g. `quant`'s
    result reduces to a 2^(w-1) or 2^(w-2) denominator whenever the value has
    trailing zero bits) -- pulling `.numerator` alone and pairing it with a
    hardcoded `2^-w` Hoon exponent then bakes in a value that's off by a power
    of two.  Always use this (never `.numerator` directly) when emitting a hex
    literal meant to be read back as `numerator / 2^w`."""
    assert (1 << w) % fr.denominator == 0
    return fr.numerator * ((1 << w) // fr.denominator)

def quant(mv): return quant_w(mv, WBITS)
def raw128(fr): return raw_w(fr, WBITS)

def val(p, n):
    """posit pattern -> Fraction, or None for NaR (mirrors posit_check._val)."""
    u = decode(p, n)
    if u[0] == 'z': return Fraction(0)
    if u[0] == 'n': return None
    _, s, e, a = u
    v = Fraction(a) * (Fraction(2) ** e)
    return -v if s else v

def nar(n): return 1 << (n - 1)

def hexn(v, n):
    return f"0x{v:0{(n+3)//4}x}"

def round_half_even(fr):
    """Fraction -> nearest int, ties to even (matches +bit's round-to-nearest-even)."""
    fl = fr.numerator // fr.denominator
    rem = fr - fl
    if rem < Fraction(1, 2): return fl
    if rem > Fraction(1, 2): return fl + 1
    return fl if fl % 2 == 0 else fl + 1

def true_pattern(true_mpf, n):
    """mpmath high-precision true value -> the CORRECTLY-ROUNDED posit pattern
    (via the same single-round `ref_value_encode` /lib/unum uses), including
    correct saturation to maxpos/minpos.  Comparing bit patterns directly
    (rather than an approximate ULP distance) avoids the pathology where a
    correctly-saturated result far from `true_mpf` looks like a huge "ULP
    error" when it's actually the right answer.

    Uses `mp.frexp` to keep 300 bits of SIGNIFICANT precision regardless of
    magnitude, not 300 bits of fixed absolute scale -- an earlier version did
    `Fraction(round(true_mpf * 2**300), 2**300)` directly, which silently
    collapses any true_mpf smaller than ~2^-300 to an exact Fraction(0)
    *before* `ref_value_encode`'s own tiny-magnitude-saturates-to-minpos
    logic ever runs, mis-grading extreme-underflow cases as "should be 0"
    when they should be minpos (posits have no underflow-to-zero). Caught by
    exhaustive posit8 sweep hitting exp(-16777216): correct output is minpos
    (0x01), the old version wrongly graded that as a 1-ULP error vs 0x00."""
    if true_mpf == 0:
        return ref_value_encode(Fraction(0), n)
    man, exp = mp.frexp(true_mpf)          # true_mpf = man * 2**exp, 0.5<=|man|<1
    SIGBITS = 300
    scaled = int(mp.nint(man * (mp.mpf(2) ** SIGBITS)))
    fr = Fraction(scaled, 1 << SIGBITS) * (Fraction(2) ** exp)
    return ref_value_encode(fr, n)

def ulp_distance(a, b, n):
    """Integer distance between two posit bit patterns' two's-complement
    ordering (posit ordering == raw two's-complement order, /lib/unum's own
    documented property) -- a robust "how many ULPs apart" for near-misses."""
    def tc(p): return p - (1 << n) if p >= (1 << (n - 1)) else p
    return abs(tc(a) - tc(b))

# ==================================== exp ====================================
# x = k*ln2 + r  (r "small"); exp(x) = 2^k * poly(r).  No hi/lo split needed
# (see module docstring); k via a stored approximate 1/ln2 (avoids an exact
# division), r = x - k*LN2 exact given LN2 is an exact (quantized) Fraction.
# Because the final combine is `encode(sign, e + k, a, n)` and `encode`
# already saturates out-of-range exponents to maxpos/minpos, overflow and
# underflow of `exp` fall out for free -- no separate guard chain needed.

LN2 = quant(mp.log(2))
INVLN2 = quant(1 / mp.log(2))

def gen_exp_coeffs(deg):
    half = mp.log(2) / 2
    cs = mp.chebyfit(lambda r: mp.e ** r, [-half, half], deg + 1)   # highest-first
    return [quant(c) for c in reversed(cs)]                          # ascending c0..cN

def horner(coeffs, r):             # ascending Fraction coeffs, exact Fraction Horner
    acc = coeffs[-1]
    for c in reversed(coeffs[:-1]):
        acc = acc * r + c
    return acc

def exp_exact(x, coeffs):
    """x: Fraction (exact decoded posit value) -> Fraction (exact algorithm
    output, before the single final `encode`)."""
    k = round_half_even(x * INVLN2)
    r = x - k * LN2
    p = horner(coeffs, r)
    return p * (Fraction(2) ** k)

def exp_ref(pattern, n, coeffs):
    x = val(pattern, n)
    if x is None: return nar(n)
    if x == 0: return encode(False, 0, 1, n)          # exp(0) = 1, exact
    return ref_value_encode(exp_exact(x, coeffs), n)

def check_exp():
    deg = 7                            # smallest degree that is correctly-rounded
                                        # (0 ULP max) at p8/p16/p32 over x in [-40,40]
    coeffs = gen_exp_coeffs(deg)
    print(f"# unum exp: k*ln2+r reduction (no hi/lo split; exact Fraction arithmetic)")
    print(f"# degree-{deg} minimax poly, WBITS={WBITS}")
    print(f"LN2    = {hexn(raw128(LN2), 8*((WBITS+7)//8+1))}  (den 2^{WBITS})")
    print(f"INVLN2 = {hexn(raw128(INVLN2), 8*((WBITS+7)//8+1))}  (den 2^{WBITS})")
    print("coeffs (ascending, c0..c%d), numerator over 2^%d:" % (deg, WBITS))
    for i, c in enumerate(coeffs):
        print(f"  c{i:<2} = {hexn(raw128(c), 8*((WBITS+7)//8+1))}   (~{float(c):.17g})")

    # Bounded-range grid (not a uniform sweep of raw bit patterns): for |x|
    # beyond ~40, exp trivially saturates to maxpos/minpos regardless of the
    # polynomial, so a raw-pattern sweep at p16/p32 wastes almost all its
    # mpmath calls on astronomically large/small posits with nothing to test.
    grid = [Fraction(t, 100) for t in range(-4000, 4001, 7)]
    for wname, n in WIDTHS:
        worst = 0; xw = None; bad = 0; tested = 0
        for xv in grid:
            p = ref_value_encode(xv, n)
            x = val(p, n)
            if x is None or x == 0: continue
            got = exp_ref(p, n, coeffs)
            true_v = mp.e ** (mp.mpf(x.numerator) / mp.mpf(x.denominator))
            want = true_pattern(true_v, n)
            tested += 1
            d = ulp_distance(got, want, n)
            if d > worst: worst, xw = d, x
            if d > 1: bad += 1         # not faithfully rounded (>1 ULP off)
        print(f"  {wname} (n={n}): max {worst} ULP at x={float(xw) if xw is not None else None}; "
              f"{bad} non-faithful (>1 ULP) of {tested} tested (x in [-40,40] grid)")

    print("edge cases (posit8):")
    for label, p in [('0', 0), ('NaR', nar(8)), ('1', encode(False,0,1,8)),
                      ('-1', encode(True,0,1,8)), ('10', ref_value_encode(Fraction(10), 8)),
                      ('50', ref_value_encode(Fraction(50), 8)), ('100', ref_value_encode(Fraction(100), 8))]:
        o = exp_ref(p, 8, coeffs)
        print(f"  exp({label:>4}) in={hexn(p,8)} -> out={hexn(o,8)}")

    print("expected bit patterns (rpb/rph/rps), for the Hoon test suite:")
    for wname, n in WIDTHS:
        for label, xv in [('0', Fraction(0)), ('1', Fraction(1)), ('-1', Fraction(-1)),
                            ('half', Fraction(1,2)), ('10', Fraction(10)), ('50', Fraction(50)),
                            ('100', Fraction(100))]:
            p = ref_value_encode(xv, n) if xv != 0 else 0
            o = exp_ref(p, n, coeffs)
            print(f"  {wname} exp({label:>4}) in={hexn(p,n)} -> out={hexn(o,n)}")
    return coeffs

# ============================ g-layer simulation ============================
# A literal Python transliteration of the Hoon `up` triple (s=non-negative?,
# e=@s, a=@u) and the shared +gmul/+gadd/+gneg/+gsub/+gdiv/+g-round helpers in
# unum.hoon -- used from here on (instead of pure-Fraction math) so the oracle
# and the Hoon are bit-for-bit the SAME algorithm, including the one place
# that isn't exact: +gdiv truncates to a fixed number of guard bits rather
# than doing unbounded-precision rational division (Hoon has no native
# rational type).  GDIV_GUARD is chosen generous enough that this truncation
# is far below any posit8/16/32 ULP -- verified empirically below, not just
# assumed.

def g_of_pattern(p, n):
    """decode() output -> (s, e, a) triple, or None for NaR/zero (caller
    handles those separately, same as the Hoon arms do before touching g)."""
    d = decode(p, n)
    if d[0] in ('z', 'n'): return d[0]
    _, neg, e, a = d
    return (not neg, e, a)              # s = non-negative, matching up.s

def gmul(x, y):
    sx, ex, ax = x; sy, ey, ay = y
    return (sx == sy, ex + ey, ax * ay)

def gadd(x, y):
    sx, ex, ax = x; sy, ey, ay = y
    if ax == 0: return y
    if ay == 0: return x
    emin = min(ex, ey)
    s1 = ax << (ex - emin); s2 = ay << (ey - emin)
    if sx == sy: return (sx, emin, s1 + s2)
    if s1 > s2: return (sx, emin, s1 - s2)
    if s2 > s1: return (sy, emin, s2 - s1)
    return (True, 0, 0)

def gneg(x):
    s, e, a = x; return (not s, e, a)

def gsub(x, y): return gadd(x, gneg(y))

def gdiv(x, y, g):
    """x/y truncated to g significant bits of EXTRA precision beyond x's own
    (matches the planned Hoon +gdiv).  The shift is `g + y.bit_length()`, not
    a bare `g` -- gmul/gadd never normalize (exact big-integer combine, no
    truncation), so a Horner-accumulated x/y pair can have huge and wildly
    mismatched raw significand bit-lengths with no relation to the actual
    magnitude of x or y (e.g. tan's sin(ax)/cos(ax), each the result of a
    degree-6 Horner chain).  A bare `(x<<g)//y` silently produces a
    ZERO-BIT (or garbage) quotient whenever x's raw bit-length happens to
    undershoot y's by more than g -- caught via tan_g giving 0 at x=11
    (correct answer -225.95) before this fix.  Scaling by y's bit-length
    guarantees g bits of real precision in the quotient regardless of either
    operand's unnormalized size."""
    sx, ex, ax = x; sy, ey, ay = y
    if ax == 0: return (True, 0, 0)
    shift = g + ay.bit_length()
    q = (ax << shift) // ay
    return (sx == sy, ex - ey - shift, q)

def g_round(x):
    """up (as a triple; caller passes (True,0,0) for zero) -> int, nearest
    integer, ties away from zero -- matches the planned Hoon +g-round."""
    s, e, a = x
    if a == 0: return 0
    if e >= 0:
        mag = a << e
    else:
        sh = -e
        q = a >> sh
        rm = a & ((1 << sh) - 1)
        half = 1 << (sh - 1)
        mag = q + 1 if rm >= half else q
    return mag if s else -mag

def g_to_value(x):                      # triple -> Fraction, for sanity prints
    s, e, a = x
    v = Fraction(a) * (Fraction(2) ** e) if e >= 0 else Fraction(a, 1 << -e)
    return v if s else -v

LN2_G = (True, -128, raw128(LN2))
ONE_G = (True, 0, 1)
TWO_G = (True, 1, 1)

def gpoly(cs, z):
    """Horner over g-triples, cs highest-degree-first (mirrors the Hoon
    ?>(^ cs) / i.cs / t.cs loop in +exp and the planned +lr)."""
    acc = cs[0]
    for c in cs[1:]:
        acc = gadd(gmul(acc, z), c)
    return acc

def exp_g(pattern, n, coeffs_g):
    """g-layer re-simulation of +exp (coeffs_g: list of g-triples, highest
    degree first) -- used to self-check the Hoon-bound hex constants BEFORE
    spending an on-ship round-trip, by construction identical arithmetic to
    what unum.hoon's +exp actually does."""
    g = g_of_pattern(pattern, n)
    if g == 'n': return nar(n)
    if g == 'z': return encode(False, 0, 1, n)
    invln2_g = (True, -128, raw128(INVLN2))
    k = g_round(gmul(g, invln2_g))
    kup = (k >= 0, 0, abs(k))
    r = gsub(g, gmul(kup, LN2_G))
    p = gpoly(coeffs_g, r)
    sp, ep, ap = p
    return encode(not sp, ep + k, ap, n)

# ==================================== log ====================================
# x = m * 2^E, m in [1,2) (free from +sea's own [a,e] split -- no subnormal
# pre-scale needed, unlike /lib/math).  f=m-1, s=f/(2+f)=f/(m+1), z=s*s;
# log(m) = 2*atanh(s) = 2*s*(1 + z*P(z)), P(z) = 1/3 + z*(1/5 + z*(1/7+...)) --
# the EXACT atanh Taylor series (no minimax fit needed: z < 1/9 here, so a
# short Taylor tail already converges far past posit32's precision).
# log(x) = E*ln2 + log(m); a genuine (truncating) division shows up in `s`,
# unlike exp/its polynomial which needed none.

GDIV_GUARD = 160                        # guard bits for the one interior division

def gen_log_coeffs(deg):
    """c_k = 1/(2k+1) for k=1..deg, as g-triples, highest-degree(=k) first."""
    return [(True, -128, raw128(quant(mp.mpf(1) / (2 * k + 1))))
            for k in range(deg, 0, -1)]

def lr_g(g):
    """g-triple (x, x>0) -> (E:int, logm:g-triple), mirrors the planned +lr."""
    s, e, a = g
    lead = a.bit_length() - 1
    mup = (True, -lead, a)                       # m = a * 2^-lead, in [1,2)
    bige = lead + e
    f = gsub(mup, ONE_G)
    s_ = gdiv(f, gadd(mup, ONE_G), GDIV_GUARD)    # s = f/(m+1)
    z = gmul(s_, s_)
    poly = gpoly(LOG_COEFFS, z)
    logm = gmul(gmul(TWO_G, s_), gadd(ONE_G, gmul(z, poly)))
    return bige, logm

def log_g(pattern, n):
    g = g_of_pattern(pattern, n)
    if g == 'n': return nar(n)
    if g == 'z': return nar(n)                    # log(0) = -inf; posits have none -> NaR
    if not g[0]: return nar(n)                    # negative -> NaR
    bige, logm = lr_g(g)
    eup = (bige >= 0, 0, abs(bige))
    result = gadd(gmul(eup, LN2_G), logm)
    sp, ep, ap = result
    return encode(not sp, ep, ap, n)

LOG_COEFFS = None                                 # set by check_log() once degree is chosen

def check_log():
    global LOG_COEFFS
    #  z < 1/9 here, so a fixed Taylor tail converges geometrically -- BUT
    #  log(x) = E*ln2 + log(m) cancels heavily whenever x is near a power of
    #  2 (both terms ~0.69 in magnitude, result can be tiny), so log(m) needs
    #  much better than posit32's ~2^-28 relative accuracy to survive that
    #  cancellation.  Empirically (dense sweep near 1/2/0.5, see git history
    #  of this file for the search): degree 16 is the smallest that's
    #  correctly-rounded (0 ULP) at p8/p16/p32; degree <=12 leaves posit32
    #  off by up to hundreds of thousands of ULP right at the cancellation
    #  point (e.g. x just below 1.0).
    deg = 16
    LOG_COEFFS = gen_log_coeffs(deg)
    print(f"# unum log: mantissa/exponent split (free from +sea) + atanh series")
    print(f"# degree-{deg} EXACT Taylor tail (no minimax fit needed), "
          f"GDIV_GUARD={GDIV_GUARD}, WBITS={WBITS}")
    for i, c in enumerate(LOG_COEFFS):
        print(f"  c_{deg-i:<2}= 1/{2*(deg-i)+1:<3} = {hexn(c[2], 8*((WBITS+7)//8+1))}")

    #  Broad grid PLUS dense sampling right around powers of 2 (1, 2, 1/2) --
    #  the catastrophic-cancellation worst case for this E*ln2+log(m) split.
    broad = [Fraction(t, 1000) for t in range(1, 200001, 101)]
    near = [Fraction(1) + Fraction(t, 1 << 24) for t in range(-2000, 2001)]
    near += [Fraction(2) + Fraction(t, 1 << 24) for t in range(-2000, 2001)]
    near += [Fraction(1, 2) + Fraction(t, 1 << 24) for t in range(-2000, 2001)]
    grid = broad + near
    for wname, n in WIDTHS:
        worst = 0; xw = None; bad = 0; tested = 0
        for xv in grid:
            p = ref_value_encode(xv, n)
            x = val(p, n)
            if x is None or x == 0: continue
            got = log_g(p, n)
            true_v = mp.log(mp.mpf(x.numerator) / mp.mpf(x.denominator))
            want = true_pattern(true_v, n)
            tested += 1
            d = ulp_distance(got, want, n)
            if d > worst: worst, xw = d, x
            if d > 1: bad += 1
        print(f"  {wname} (n={n}): max {worst} ULP at x={float(xw) if xw is not None else None}; "
              f"{bad} non-faithful (>1 ULP) of {tested} tested (broad + near-power-of-2 grid)")

    print("edge/domain cases (posit8):")
    for label, p in [('0', 0), ('NaR', nar(8)), ('1', encode(False,0,1,8)),
                      ('-1', encode(True,0,1,8)), ('2', ref_value_encode(Fraction(2), 8)),
                      ('100', ref_value_encode(Fraction(100), 8))]:
        o = log_g(p, 8)
        print(f"  log({label:>4}) in={hexn(p,8)} -> out={hexn(o,8)}")

    print("expected bit patterns (rpb/rph/rps), for the Hoon test suite:")
    for wname, n in WIDTHS:
        for label, xv in [('1', Fraction(1)), ('2', Fraction(2)), ('half', Fraction(1,2)),
                            ('10', Fraction(10)), ('100', Fraction(100))]:
            p = ref_value_encode(xv, n)
            o = log_g(p, n)
            print(f"  {wname} log({label:>4}) in={hexn(p,n)} -> out={hexn(o,n)}")
    return LOG_COEFFS

# ============================= log-2 / log-10 ================================
# log_b(x) = E*log_b(2) + log(m)/ln(b), combined directly from +lr's [E, logm]
# -- NOT `log(x)/log_b` (that would double-round through a posit-rounded
# log2/log10 constant).  Needs LOG_COEFFS set (see check_log/check_log2log10).

LOG10_2_G = (True, -128, raw128(quant(mp.log(2) / mp.log(10))))
INVLN10_G = (True, -128, raw128(quant(1 / mp.log(10))))

def log2_g(pattern, n):
    g = g_of_pattern(pattern, n)
    if g in ('n', 'z') or not g[0]: return nar(n)
    bige, logm = lr_g(g)
    eup = (bige >= 0, 0, abs(bige))
    invln2_g = (True, -128, raw128(INVLN2))
    sp, ep, ap = gadd(eup, gmul(logm, invln2_g))
    return encode(not sp, ep, ap, n)

def log10_g(pattern, n):
    g = g_of_pattern(pattern, n)
    if g in ('n', 'z') or not g[0]: return nar(n)
    bige, logm = lr_g(g)
    eup = (bige >= 0, 0, abs(bige))
    sp, ep, ap = gadd(gmul(eup, LOG10_2_G), gmul(logm, INVLN10_G))
    return encode(not sp, ep, ap, n)

def check_log2log10():
    global LOG_COEFFS
    LOG_COEFFS = gen_log_coeffs(16)
    print(f"log10-2-wide = {hexn(raw128(quant(mp.log(2)/mp.log(10))), 8*((WBITS+7)//8+1))}")
    print(f"invln10-wide = {hexn(raw128(quant(1/mp.log(10))), 8*((WBITS+7)//8+1))}")
    #  NOTE: this function previously only printed bit patterns for a handful
    #  of round values with no ULP sweep at all -- added below (mirrors
    #  check_log's own broad + near-power-of-2 grid) since log-2/log-10 were
    #  otherwise the only unum transcendentals never actually measured.
    broad = [Fraction(t, 1000) for t in range(1, 200001, 101)]
    near = [Fraction(1) + Fraction(t, 1 << 24) for t in range(-2000, 2001)]
    near += [Fraction(2) + Fraction(t, 1 << 24) for t in range(-2000, 2001)]
    near += [Fraction(1, 2) + Fraction(t, 1 << 24) for t in range(-2000, 2001)]
    grid = broad + near
    for name, fn, truefn in [('log2', log2_g, lambda xm: mp.log(xm) / mp.log(2)),
                              ('log10', log10_g, lambda xm: mp.log(xm) / mp.log(10))]:
        for wname, n in WIDTHS:
            worst = 0; xw = None; bad = 0; tested = 0
            for xv in grid:
                p = ref_value_encode(xv, n)
                x = val(p, n)
                if x is None or x == 0: continue
                got = fn(p, n)
                true_v = truefn(mp.mpf(x.numerator) / mp.mpf(x.denominator))
                want = true_pattern(true_v, n)
                tested += 1
                d = ulp_distance(got, want, n)
                if d > worst: worst, xw = d, x
                if d > 1: bad += 1
            print(f"  {name} {wname} (n={n}): max {worst} ULP at x={float(xw) if xw is not None else None}; "
                  f"{bad} non-faithful (>1 ULP) of {tested} tested (broad + near-power-of-2 grid)")
    for name, fn in [('log2', log2_g), ('log10', log10_g)]:
        for wname, n in WIDTHS:
            for label, xv in [('1', Fraction(1)), ('2', Fraction(2)), ('half', Fraction(1, 2)),
                                ('10', Fraction(10)), ('100', Fraction(100))]:
                p = ref_value_encode(xv, n)
                o = fn(p, n)
                print(f"  {wname} {name}({label:>4}) in={hexn(p,n)} -> out={hexn(o,n)}")

# =================================== sin/cos ===================================
# Quarter-turn reduction on |x| (sin is odd, cos is even -- reduce sign
# separately so the reduced q=round(|x|*2/pi) is always >=0, no negative-mod
# bookkeeping needed).  r = |x| - q*(pi/2), single wide pi/2 constant (no
# multi-part split -- see module docstring).  Kernels are the EXACT Taylor
# series in z=r*r (r in [-pi/4,pi/4]), same "no minimax needed, just enough
# terms" approach as +log's atanh series.
#
# UNLIKE exp/log, large-argument reduction error scales with the ARGUMENT
# (error in pi/2 gets multiplied by q ~ x*2/pi), so WBITS=128 isn't enough to
# stay correctly-rounded near posit32's maxpos (~2^120) -- /lib/math hits the
# same wall on hardware floats and documents a bounded "faithful" range
# (e.g. "|x|<~500" for @rs).  Because Hoon atoms are arbitrary-precision this
# is just a bigger constant, not new machinery: TRIG_WBITS=200 gives enough
# margin (200 - 120 = 80 spare bits) to stay correctly-rounded across
# posit32's ENTIRE dynamic range (verified below, exponential sweep to 2^120)
# -- no documented limitation needed here, unlike /lib/math.

TRIG_WBITS = 200
PI2 = quant_w(mp.pi / 2, TRIG_WBITS)
INVPI2 = quant_w(2 / mp.pi, TRIG_WBITS)

def gen_sincos_coeffs(deg):
    """(sin_coeffs, cos_coeffs), each highest-degree-first g-triples.
    sin(r) = r*(1 + z*Psin(z)), Psin c_k=(-1)^k/(2k+1)! for k=1..deg.
    cos(r) = 1 + z*Pcos(z),      Pcos c_k=(-1)^k/(2k)!   for k=1..deg."""
    sin_cs = [(True if k % 2 == 0 else False, -128,
               raw128(quant(mp.mpf(1) / mp.factorial(2 * k + 1))))
              for k in range(deg, 0, -1)]
    cos_cs = [(True if k % 2 == 0 else False, -128,
               raw128(quant(mp.mpf(1) / mp.factorial(2 * k))))
              for k in range(deg, 0, -1)]
    return sin_cs, cos_cs

SIN_COEFFS = COS_COEFFS = None            # set by check_trig()
PI2_G = None; INVPI2_G = None; ONE_G = (True, 0, 1)

def sc_g(g):
    """g-triple ax (>=0, %p only -- caller handles NaR/zero) -> (sin(ax),
    cos(ax)) as UNROUNDED g-triples."""
    q = g_round(gmul(g, INVPI2_G))
    qn = abs(q)                                       # q always >=0 here (ax>=0, invpi2>0)
    r = gsub(g, gmul((True, 0, qn), PI2_G))
    z = gmul(r, r)
    sink = gmul(r, gadd(ONE_G, gmul(z, gpoly(SIN_COEFFS, z))))
    cosk = gadd(ONE_G, gmul(z, gpoly(COS_COEFFS, z)))
    m = qn % 4
    if m == 0: return sink, cosk
    if m == 1: return cosk, gneg(sink)
    if m == 2: return gneg(sink), gneg(cosk)
    return gneg(cosk), sink

def sin_g(pattern, n):
    g = g_of_pattern(pattern, n)
    if g == 'n': return nar(n)
    if g == 'z': return 0
    neg_in = not g[0]                                  # was x negative?
    ax = (True, g[1], g[2])
    sn, _ = sc_g(ax)
    sp, ep, ap = sn
    final_sp = (not sp) if neg_in else sp               # sin is odd: flip if x was negative
    return encode(not final_sp, ep, ap, n)

TAN_GUARD = 160

def tan_g(pattern, n):
    g = g_of_pattern(pattern, n)
    if g == 'n': return nar(n)
    if g == 'z': return 0
    neg_in = not g[0]
    ax = (True, g[1], g[2])
    sn, cs = sc_g(ax)
    if cs[2] == 0: return nar(n)                        # cos(ax) exactly 0 (measure-zero, but guard anyway)
    raw = gdiv(sn, cs, TAN_GUARD)
    sp, ep, ap = raw
    final_sp = (not sp) if neg_in else sp                # tan is odd, same sign rule as sin
    return encode(not final_sp, ep, ap, n)

def cos_g(pattern, n):
    g = g_of_pattern(pattern, n)
    if g == 'n': return nar(n)
    if g == 'z': return encode(False, 0, 1, n)          # cos(0) = 1
    ax = (True, g[1], g[2])                             # cos is even: sign doesn't matter
    _, cs = sc_g(ax)
    sp, ep, ap = cs
    return encode(not sp, ep, ap, n)

def check_trig():
    global SIN_COEFFS, COS_COEFFS, PI2_G, INVPI2_G
    deg = 6                      # smallest degree that's correctly-rounded at
                                  # p8/16/32; z up to (pi/4)^2 ~ 0.617
    SIN_COEFFS, COS_COEFFS = gen_sincos_coeffs(deg)
    PI2_G = (True, -TRIG_WBITS, raw_w(PI2, TRIG_WBITS))
    INVPI2_G = (True, -TRIG_WBITS, raw_w(INVPI2, TRIG_WBITS))
    print(f"# unum sin/cos: quarter-turn |x| reduction + degree-{deg} EXACT Taylor kernels")
    print(f"# TRIG_WBITS={TRIG_WBITS} (wider than WBITS={WBITS} -- see module docstring: "
          f"large-x reduction error scales with the argument, not just the constant's own precision)")
    print(f"PI2    = {hexn(raw_w(PI2, TRIG_WBITS), 52)}")
    print(f"INVPI2 = {hexn(raw_w(INVPI2, TRIG_WBITS), 52)}")

    broad = [Fraction(t, 1000) for t in range(0, 200001, 101)]
    # dense near quadrant boundaries (k*pi/4) -- cancellation-ish worst case
    near = []
    for k in range(0, 9):
        c = quant(mp.pi * k / 4)
        near += [c + Fraction(t, 1 << 40) for t in range(-500, 501)]
    grid = broad + near
    #  NOTE: tan_g is defined (used internally, e.g. the GDIV_GUARD tuning
    #  comment above) but had NO dedicated ULP check anywhere in this file --
    #  added here (same grid as sin/cos) since it's otherwise the only unum
    #  trig arm never actually measured.
    for fname, fn, truefn in [('sin', sin_g, mp.sin), ('cos', cos_g, mp.cos), ('tan', tan_g, mp.tan)]:
        for wname, n in WIDTHS:
            worst = 0; xw = None; bad = 0; tested = 0
            for xv in grid:
                for sgn in (1, -1):
                    p = ref_value_encode(sgn * xv, n)
                    x = val(p, n)
                    if x is None: continue
                    got = fn(p, n)
                    true_v = truefn(mp.mpf(x.numerator) / mp.mpf(x.denominator))
                    if fname == 'tan' and (not mp.isfinite(true_v) or abs(true_v) > mp.mpf('1e15')): continue
                    want = true_pattern(true_v, n)
                    tested += 1
                    d = ulp_distance(got, want, n)
                    if d > worst: worst, xw = d, x
                    if d > 1: bad += 1
            print(f"  {fname} {wname} (n={n}): max {worst} ULP at x={float(xw) if xw is not None else None}; "
                  f"{bad} non-faithful (>1 ULP) of {tested} tested")

    #  Exponential-magnitude sweep to posit32's maxpos (~2^120) -- confirms
    #  TRIG_WBITS=200 needs no documented large-x limitation (unlike math.hoon).
    import random
    random.seed(5)
    mp.mp.dps = 120
    worst_huge = 0; huge_x = None
    for mag_exp in range(0, 121, 4):
        base = mp.mpf(2) ** mag_exp
        for _ in range(20):
            xf = base * (1 + random.random())
            xv = Fraction(int(xf * 10**6), 10**6)
            for sgn in (1, -1):
                p = ref_value_encode(sgn * xv, 32)
                x = val(p, 32)
                if x is None: continue
                for fn, truefn in [(sin_g, mp.sin), (cos_g, mp.cos)]:
                    got = fn(p, 32)
                    tv = truefn(mp.mpf(x.numerator) / mp.mpf(x.denominator))
                    want = true_pattern(tv, 32)
                    d = ulp_distance(got, want, 32)
                    if d > worst_huge: worst_huge, huge_x = d, x
    mp.mp.dps = 60
    print(f"  huge-x sweep (posit32, |x| up to ~2^120): max {worst_huge} ULP"
          f" at x~{float(huge_x) if huge_x else None}")

    print("expected bit patterns (rpb/rph/rps), for the Hoon test suite:")
    for wname, n in WIDTHS:
        for label, xv in [('0', Fraction(0)), ('1', Fraction(1)), ('-1', Fraction(-1)),
                            ('pi2', quant(mp.pi/2)), ('pi', quant(mp.pi)), ('10', Fraction(10))]:
            p = ref_value_encode(xv, n) if xv != 0 else 0
            print(f"  {wname} sin({label:>4}) in={hexn(p,n)} -> out={hexn(sin_g(p,n),n)}"
                  f"   cos({label:>4}) -> out={hexn(cos_g(p,n),n)}")

# ==================================== atan ====================================
# fdlibm breakpoint reduction (7/16, 11/16, 19/16, 39/16 -- exact rational
# comparisons, no quantization needed for the thresholds themselves): pick a
# breakpoint bp in {0, atan(1/2), pi/4, atan(3/2), pi/2} via tan-subtraction
# xr = (ax-bp_val)/(1+ax*bp_val), so atan(ax) = bp_angle + atan(xr) with |xr|
# always small.  atan(xr) itself via the EXACT odd Taylor series in z=xr*xr
# (no minimax fit needed, same "no hi/lo split, no argument-scaled error"
# situation as +log -- the breakpoint angle is a fixed additive constant, not
# multiplied by a growing integer the way sin/cos's q*pi2 is).

BP_VALS = [Fraction(0), Fraction(1, 2), Fraction(1), Fraction(3, 2)]     # tan(breakpoint angle)
ATAN_COEFFS = None                                                        # set by check_atan()

def atan_core(ax):
    """ax: g-triple, x>=0, %p only -> raw unrounded atan(ax) as a g-triple."""
    axv = g_to_value(ax)
    if axv < Fraction(7, 16):
        band, bp_val = 0, Fraction(0)
    elif axv < Fraction(11, 16):
        band, bp_val = 1, BP_VALS[1]
    elif axv < Fraction(19, 16):
        band, bp_val = 2, BP_VALS[2]
    elif axv < Fraction(39, 16):
        band, bp_val = 3, BP_VALS[3]
    else:
        band = 4; bp_val = None
    if band == 0:
        xr = ax
    elif band == 4:
        xr = gneg(gdiv(ONE_G, ax, 160))
    else:
        bp_g = (True, -128, raw128(quant(mp.mpf(bp_val.numerator) / bp_val.denominator)))
        num = gsub(ax, bp_g)
        den = gadd(ONE_G, gmul(ax, bp_g))
        xr = gdiv(num, den, 160)
    z = gmul(xr, xr)
    series = gadd(ONE_G, gmul(z, gpoly(ATAN_COEFFS, z)))
    at_xr = gmul(xr, series)
    bp_angle_g = BP_ANGLES[band]
    return gadd(bp_angle_g, at_xr)

def gen_atan_coeffs(deg):
    """c_k = (-1)^k/(2k+1) for k=1..deg, highest-degree first."""
    return [(k % 2 == 0, -128, raw128(quant(mp.mpf(1) / (2 * k + 1))))
            for k in range(deg, 0, -1)]

def atan_g(pattern, n):
    g = g_of_pattern(pattern, n)
    if g == 'n': return nar(n)
    if g == 'z': return 0
    neg_in = not g[0]
    ax = (True, g[1], g[2])
    raw = atan_core(ax)
    sp, ep, ap = raw
    final_sp = (not sp) if neg_in else sp
    return encode(not final_sp, ep, ap, n)

def check_atan():
    global ATAN_COEFFS, BP_ANGLES
    deg = 13                     # smallest degree that's correctly-rounded at p8/16/32
    ATAN_COEFFS = gen_atan_coeffs(deg)
    BP_ANGLES = [(True, 0, 0),
                 (True, -128, raw128(quant(mp.atan(mp.mpf(1) / 2)))),
                 (True, -128, raw128(quant(mp.pi / 4))),
                 (True, -128, raw128(quant(mp.atan(mp.mpf(3) / 2)))),
                 (True, -128, raw128(quant(mp.pi / 2)))]
    print(f"# unum atan: fdlibm breakpoint reduction + degree-{deg} EXACT Taylor kernel")

    broad = [Fraction(t, 1000) for t in range(0, 200001, 101)]
    near = []
    for bpv in [Fraction(7, 16), Fraction(11, 16), Fraction(19, 16), Fraction(39, 16)]:
        near += [bpv + Fraction(t, 1 << 30) for t in range(-500, 501)]
    grid = broad + near
    for wname, n in WIDTHS:
        worst = 0; xw = None; bad = 0; tested = 0
        for xv in grid:
            for sgn in (1, -1):
                p = ref_value_encode(sgn * xv, n)
                x = val(p, n)
                if x is None: continue
                got = atan_g(p, n)
                true_v = mp.atan(mp.mpf(x.numerator) / mp.mpf(x.denominator))
                want = true_pattern(true_v, n)
                tested += 1
                d = ulp_distance(got, want, n)
                if d > worst: worst, xw = d, x
                if d > 1: bad += 1
        print(f"  atan {wname} (n={n}): max {worst} ULP at x={float(xw) if xw is not None else None}; "
              f"{bad} non-faithful (>1 ULP) of {tested} tested")

    print("expected bit patterns (rpb/rph/rps), for the Hoon test suite:")
    for wname, n in WIDTHS:
        for label, xv in [('0', Fraction(0)), ('1', Fraction(1)), ('-1', Fraction(-1)),
                            ('half', Fraction(1, 2)), ('10', Fraction(10))]:
            p = ref_value_encode(xv, n) if xv != 0 else 0
            print(f"  {wname} atan({label:>4}) in={hexn(p,n)} -> out={hexn(atan_g(p,n),n)}")

# ================================= asin / acos =================================
# Compose from already-correctly-rounded posit ops (+mul/+sub/+sqt/+div, all
# existing and unchanged) plus the new precise +atan, mirroring /lib/unum's
# EXISTING identities (asin(x)=atan(x/sqrt(1-x^2)), acos(x)=atan(sqrt(1-x^2)/x)
# for x>0, pi-acos(-x) for x<0) rather than a dedicated rational P/Q kernel
# like /lib/math's @rd/@rs -- much simpler, and the naive version already only
# claimed ~1 ULP here (AGM+sqrt rounding).  This composition is FAITHFUL
# (<=7 ULP observed, see check_ainv() below) but not always correctly rounded
# -- the two intermediate roundings (1-x^2, sqrt, the ratio) each cost a
# fraction of a ULP, compounding near |x|~1 where asin/acos's derivative
# blows up.  Accepted as consistent with /lib/math's own ~0.8-1 ULP bound on
# asin/acos; a dedicated kernel would remove this but wasn't judged worth the
# extra implementation complexity for this pass.

def _one_p(n): return ref_value_encode(Fraction(1), n)
def _negone_p(n): return ref_value_encode(Fraction(-1), n)

def asin_g(pat, n):
    d = decode(pat, n)
    if d[0] == 'n': return nar(n)
    if pat == _one_p(n): return true_pattern(mp.pi / 2, n)
    if pat == _negone_p(n): return true_pattern(-mp.pi / 2, n)
    if d[0] == 'z': return 0
    x = val(pat, n)
    if x is None or abs(x) > 1: return nar(n)
    t = pc_sub(_one_p(n), pc_mul(pat, pat, n), n)
    s = my_sqrt(t, n)
    q = pc_div(pat, s, n)
    return atan_g(q, n)

def acos_g(pat, n):
    d = decode(pat, n)
    if d[0] == 'n': return nar(n)
    if pat == _one_p(n): return 0
    if pat == _negone_p(n): return true_pattern(mp.pi, n)
    if d[0] == 'z': return true_pattern(mp.pi / 2, n)
    x = val(pat, n)
    if x is None or abs(x) > 1: return nar(n)
    t = pc_sub(_one_p(n), pc_mul(pat, pat, n), n)
    s = my_sqrt(t, n)
    ax = pat if x > 0 else pc_neg(pat, n)
    q = pc_div(s, ax, n)
    at = atan_g(q, n)
    if x > 0: return at
    return pc_sub(true_pattern(mp.pi, n), at, n)      # acos(x) = pi - acos(-x)

def check_ainv():
    global ATAN_COEFFS, BP_ANGLES
    ATAN_COEFFS = gen_atan_coeffs(13)
    BP_ANGLES = [(True, 0, 0),
                 (True, -128, raw128(quant(mp.atan(mp.mpf(1) / 2)))),
                 (True, -128, raw128(quant(mp.pi / 4))),
                 (True, -128, raw128(quant(mp.atan(mp.mpf(3) / 2)))),
                 (True, -128, raw128(quant(mp.pi / 2)))]
    broad = [Fraction(t, 1000) for t in range(-999, 1000, 7)]
    for name, fn, truefn in [('asin', asin_g, mp.asin), ('acos', acos_g, mp.acos)]:
        for wname, n in WIDTHS:
            worst = 0; xw = None; bad = 0; tested = 0
            for xv in broad:
                p = ref_value_encode(xv, n)
                x = val(p, n)
                if x is None or abs(x) > 1: continue
                got = fn(p, n)
                true_v = truefn(mp.mpf(x.numerator) / mp.mpf(x.denominator))
                want = true_pattern(true_v, n)
                tested += 1
                d = ulp_distance(got, want, n)
                if d > worst: worst, xw = d, x
                if d > 1: bad += 1
            print(f"  {name} {wname} (n={n}): max {worst} ULP at x={float(xw) if xw is not None else None}; "
                  f"{bad} non-faithful (>1 ULP) of {tested} tested")
    print("expected bit patterns (rpb/rph/rps), for the Hoon test suite:")
    for wname, n in WIDTHS:
        for label, xv in [('0', Fraction(0)), ('1', Fraction(1)), ('-1', Fraction(-1)),
                            ('half', Fraction(1, 2)), ('-half', Fraction(-1, 2))]:
            p = ref_value_encode(xv, n) if xv != 0 else 0
            print(f"  {wname} asin({label:>5}) in={hexn(p,n)} -> out={hexn(asin_g(p,n),n)}"
                  f"   acos({label:>5}) -> out={hexn(acos_g(p,n),n)}")

# ================================= cbrt / pow / pow-n =================================
#  All three are literal compositions of already-checked primitives in
#  lib/unum.hoon (not separate polynomial kernels):
#    +pow-n:  [x=@ p=@u] -> @   repeated `mul` (posit-rounded multiply), p times
#    +pow:    [x=@ y=@] -> @   (exp (mul y (log x)))  -- ALWAYS via exp/log, no
#                               integer-exponent fast path (unlike /lib/math's
#                               +pow, which falls through to +pow-n for y>0 int)
#    +cbrt:   @ -> @            (pow x (div one (sun 3)))  -- note this divides
#                               to get a POSIT-ROUNDED 1/3 (not an exact
#                               constant like /lib/math's), an extra rounding
#                               step /lib/math's cbrt doesn't have.
#  This section re-executes those exact compositions using the already-
#  validated exp_ref/log_g/pc_mul/pc_div (full simulation, not a derived
#  bound), against the same true_pattern single-round oracle used above.
POW_EXP_COEFFS = None                  # set by check_pow_pown_cbrt()

def pow_n_ref(x_pat, p, n):
    """mirrors +pow-n: [x=@ p=@u] -> @ (repeated posit-rounded multiply)."""
    if x_pat == nar(n): return nar(n)
    res = _one_p(n)
    for _ in range(p):
        res = pc_mul(res, x_pat, n)
    return res

def pow_ref(x_pat, y_pat, n, exp_coeffs):
    """mirrors +pow: exp(mul(y, log(x))) -- no integer-exponent fast path."""
    lx = log_g(x_pat, n)
    if lx == nar(n): return nar(n)
    my = pc_mul(y_pat, lx, n)
    return exp_ref(my, n, exp_coeffs)

def cbrt_ref(x_pat, n, exp_coeffs):
    """mirrors +cbrt: nar/0 short-circuit, x<0 -> nar, else pow(x, div(one,sun(3)))."""
    if x_pat == nar(n): return nar(n)
    if x_pat == 0: return 0
    xv = val(x_pat, n)
    if xv is not None and xv < 0: return nar(n)
    third = pc_div(_one_p(n), ref_value_encode(Fraction(3), n), n)
    return pow_ref(x_pat, third, n, exp_coeffs)

def check_pow_n():
    print("# unum pow-n: repeated posit-rounded multiply -- EXACT rational ground truth (no mpmath needed)")
    for wname, n in WIDTHS:
        for p in (2, 3, 5, 7):
            worst = 0; xw = None; bad = 0; tested = 0
            for t in range(-300, 301):
                if t == 0: continue
                xv = Fraction(t, 30)
                xp = ref_value_encode(xv, n)
                x = val(xp, n)
                if x is None: continue
                got = pow_n_ref(xp, p, n)
                want = ref_value_encode(x ** p, n)     # exact rational power -> exact oracle
                tested += 1
                d = ulp_distance(got, want, n)
                if d > worst: worst, xw = d, x
                if d > 1: bad += 1
            print(f"  pow-n(x,{p}) {wname} (n={n}): max {worst} ULP at x={float(xw) if xw is not None else None}; "
                  f"{bad} non-faithful (>1 ULP) of {tested} tested")

def check_pow():
    global POW_EXP_COEFFS, LOG_COEFFS
    LOG_COEFFS = gen_log_coeffs(16)
    POW_EXP_COEFFS = gen_exp_coeffs(7)
    print("# unum pow: exp(mul(y, log(x))) -- composed; full simulation (mpmath oracle)")
    grid_x = [Fraction(t, 10) for t in range(1, 401, 7)]
    grid_y = [Fraction(t, 10) for t in range(-100, 101, 7)]
    for wname, n in WIDTHS:
        worst = 0; xw = None; bad = 0; tested = 0
        for xv in grid_x:
            xp = ref_value_encode(xv, n)
            x = val(xp, n)
            if x is None or x <= 0: continue
            for yv in grid_y:
                yp = ref_value_encode(yv, n)
                y = val(yp, n)
                if y is None: continue
                got = pow_ref(xp, yp, n, POW_EXP_COEFFS)
                if got == nar(n): continue
                true_v = mp.mpf(x.numerator) / mp.mpf(x.denominator)
                true_v = true_v ** (mp.mpf(y.numerator) / mp.mpf(y.denominator))
                if not mp.isfinite(true_v) or true_v == 0: continue
                want = true_pattern(true_v, n)
                tested += 1
                d = ulp_distance(got, want, n)
                if d > worst: worst, xw = d, (x, y)
                if d > 1: bad += 1
        print(f"  pow {wname} (n={n}): max {worst} ULP at (x,y)={xw}; "
              f"{bad} non-faithful (>1 ULP) of {tested} tested (x in (0,40], y in [-10,10])")

def check_cbrt():
    global POW_EXP_COEFFS, LOG_COEFFS
    LOG_COEFFS = gen_log_coeffs(16)
    POW_EXP_COEFFS = gen_exp_coeffs(7)
    print("# unum cbrt: pow(x, div(one,sun(3))) -- posit-ROUNDED 1/3, not exact; full simulation (mpmath oracle)")
    print("# NOTE: unlike /lib/math's cbrt (sign(x)*exp(log|x|/3), defined for all reals), /lib/unum's")
    print("#   +cbrt explicitly returns NaR for x<0 (no sign-extraction trick) -- verified below, and the")
    print("#   ULP sweep is restricted to x>0 accordingly (comparing NaR against a negative real cube root")
    print("#   oracle would score a domain restriction as a bogus multi-billion-ULP 'error').")
    xn = ref_value_encode(Fraction(-1, 8), 8)
    print(f"  domain check: cbrt(-1/8) -> {hexn(cbrt_ref(xn, 8, POW_EXP_COEFFS), 8)}  (nar = {hexn(nar(8), 8)})")
    grid = [Fraction(t, 100) for t in range(1, 4001, 7)]
    for wname, n in WIDTHS:
        worst = 0; xw = None; bad = 0; tested = 0
        for xv in grid:
            xp = ref_value_encode(xv, n)
            x = val(xp, n)
            if x is None or x <= 0: continue
            got = cbrt_ref(xp, n, POW_EXP_COEFFS)
            true_v = mp.cbrt(mp.mpf(x.numerator) / mp.mpf(x.denominator))
            want = true_pattern(true_v, n)
            tested += 1
            d = ulp_distance(got, want, n)
            if d > worst: worst, xw = d, x
            if d > 1: bad += 1
        print(f"  cbrt {wname} (n={n}): max {worst} ULP at x={float(xw) if xw is not None else None}; "
              f"{bad} non-faithful (>1 ULP) of {tested} tested (x>0 only)")

if __name__ == '__main__':
    fn = sys.argv[1] if len(sys.argv) > 1 else 'exp'
    {'exp': check_exp, 'log': check_log, 'log2log10': check_log2log10,
     'trig': check_trig, 'atan': check_atan, 'ainv': check_ainv,
     'pow-n': check_pow_n, 'pow': check_pow, 'cbrt': check_cbrt}[fn]()
