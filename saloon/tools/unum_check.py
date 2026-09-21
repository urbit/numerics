#!/usr/bin/env python3
"""Oracle for saloon/desk/tests/lib/saloon-unum.hoon.

Recomputes every expected posit32 value in that file, correctly rounded, and
checks it against the literal the test asserts.  Run with no arguments; exits
nonzero on any mismatch.

The reference is libmath/tools/posit_check.py's exact-rational encoder
(2022 Posit Standard, es=2), fed 80-digit mpmath values -- far more precision
than a posit32 rounding decision needs.  This is independent of /lib/unum's
own implementation, so it catches the case that motivated it: the test values
were once the OLD series outputs, one ulp off for exp/cos/tan/log, and stayed
that way after numerics #71 made /lib/unum correctly rounded.

    pip install mpmath
    python3 saloon/tools/unum_check.py
"""
import os
import re
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "libmath", "tools"))
import posit_check as pc  # noqa: E402
from mpmath import mp, mpf, exp, sin, cos, tan, log, cbrt, factorial  # noqa: E402

mp.dps = 80
TEST = os.path.join(HERE, "..", "desk", "tests", "lib", "saloon-unum.hoon")


def exact(x):
    """The exact Fraction of an mpmath value (sign, mantissa, exponent)."""
    sign, man, ex, _ = mpf(x)._mpf_
    v = Fraction(int(man)) * (Fraction(2) ** int(ex))
    return -v if sign else v


half = mpf(1) / 2
#  arm name -> the exact real value the test's input produces
REFS = {
    "exp": exp(half),        # (ur 5 0x3800.0000) is posit32 0.5
    "sin": sin(half),
    "cos": cos(half),
    "tan": tan(half),
    "log": log(mpf(2)),      # (sun:rps:unum 2)
    "cbrt": cbrt(mpf(1)),    # (sun:rps:unum 1)
    "fact": factorial(3),    # (sun:rps:unum 3)
    "pown": mpf(2) ** 3,     # 2 ^ 3
}

#  sanity: the input 0x3800.0000 really is 0.5 at posit32
assert pc.ref_value_encode(Fraction(1, 2), 32) == 0x38000000

src = open(TEST).read()
found = dict(re.findall(r"\+\+  test-unum-([a-z]+)-rps\s+\(expect-eq !>\(`@`(0x[0-9a-f.]+)\)", src))
fails = 0
print(f"{'arm':6s} {'correctly rounded':>18s} {'test asserts':>14s}")
for arm, ref in REFS.items():
    want = pc.ref_value_encode(exact(ref), 32)
    if arm not in found:
        print(f"{arm:6s} {want:#018x} {'(missing)':>14s}  FAIL")
        fails += 1
        continue
    got = int(found[arm].replace(".", ""), 16)
    ok = got == want
    fails += not ok
    print(f"{arm:6s} {want:#018x} {got:#014x}  {'ok' if ok else 'FAIL'}")
extra = set(found) - set(REFS)
if extra:
    print(f"untested arms in the file: {sorted(extra)}")
    fails += len(extra)
print(f"\n{len(REFS)} checks, {fails} failures")
sys.exit(1 if fails else 0)
