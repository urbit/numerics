#!/usr/bin/env python3
"""Oracle for librand's /lib/unumrand +posit-unit (rand-spec.md section 12.4).

+posit-unit draws k raw bits u and encodes the dyadic u/2**k (the g-layer
value [%p %.y -k u]) via /lib/unum's +bit (round-to-nearest-even,
saturating) to get a uniform value on [0,1) at posit precision. This script
is the independent oracle: it reuses libmath/tools/posit_check.py's own
decode/encode reference model (already cross-checked against SoftPosit) to
derive, EXACTLY -- via Fraction arithmetic, no floating point anywhere --
the true probability every reachable posit pattern must receive under this
construction. That table is then chi-squared against empirical counts drawn
on-ship at a fixed seed.

decode()/encode() below are a verbatim copy of libmath/tools/posit_check.py's
reference model (not imported, to keep this script standalone and runnable
from any directory) -- see that file's own header for its SoftPosit
cross-check history. Everything past that point (the expected-probability
derivation, binning, chi-square) is new, specific to +posit-unit.

Usage:
  # Inspect the exact expected-probability table:
  python3 posit_unit_check.py table posit8
  python3 posit_unit_check.py table posit16

  # Chi-square ship-drawn raw patterns (one integer per line, hex 0x.. or
  # decimal) against the exact expected table, binned into equal-mass groups:
  python3 posit_unit_check.py chi2 posit8 counts.txt --bins 32
  python3 posit_unit_check.py chi2 posit16 counts.txt --bins 64

Optional: scipy, for a p-value alongside the raw statistic. Falls back to
reporting just the statistic and degrees of freedom if scipy isn't
installed.
"""
import sys
import argparse
from fractions import Fraction

try:
    from scipy.stats import chi2 as _scipy_chi2
except ImportError:
    _scipy_chi2 = None

WIDTHS = {'posit8': 8, 'posit16': 16, 'posit32': 32}
# k = 4n, per rand-spec.md 12.4: gives a constant 7-bit safety margin over
# the finest rounding-cell width in [0,1) at any n (that width sits at
# exponent -(4(n-2)+1), so k=4n exceeds it by exactly 4n - (4(n-2)+1) = 7).
KBITS = {'posit8': 32, 'posit16': 64, 'posit32': 128}

# ---- verbatim reference model (see module docstring) ----

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

def value(p, n):
    """Exact Fraction value of pattern p at width n, or None for NaR."""
    d = decode(p, n)
    if d[0] == 'z': return Fraction(0)
    if d[0] == 'n': return None
    _, neg, e, a = d
    v = Fraction(a) * (Fraction(2) ** e)
    return -v if neg else v

# ---- exact expected-probability table for +posit-unit ----

def reachable_patterns(n):
    """Nonneg-valued patterns (0 plus positive posits) with value <= 1,
    sorted by value ascending. +posit-unit can only ever produce one of
    these: the input is always in [0, 1 - 2**-k], and RNE rounding of a
    nonnegative value never flips sign or overshoots past 1.0's own
    pattern."""
    entries = []
    for p in range(1 << n):
        d = decode(p, n)
        if d[0] == 'n':
            continue
        if d[0] == 'z':
            entries.append((Fraction(0), p))
            continue
        _, neg, e, a = d
        if neg:
            continue
        v = Fraction(a) * (Fraction(2) ** e)
        if v <= 1:
            entries.append((v, p))
    entries.sort()
    return entries

def expected_probabilities(width_name):
    """Exact Fraction probability for every reachable pattern.

    For each pattern, finds the exact range of numerators u in [0, 2**k)
    that round to it, via binary search against the TRUSTED `encode`
    (the same round-to-nearest-even logic /lib/unum's +bit implements) --
    not by re-deriving midpoint/tie-break math by hand. Nonneg posit
    patterns are monotonic in both raw-integer and decoded-value order, so
    comparing raw pattern integers directly (no re-decode) is enough.
    """
    n = WIDTHS[width_name]
    k = KBITS[width_name]
    two_k = 1 << k
    entries = reachable_patterns(n)
    patterns = [p for _, p in entries]
    m = len(patterns)

    def first_u_reaching(target_pattern):
        lo, hi = 0, two_k
        while lo < hi:
            mid = (lo + hi) // 2
            out = encode(False, -k, mid, n)
            if out >= target_pattern:
                hi = mid
            else:
                lo = mid + 1
        return lo

    thresholds = [0] * (m + 1)
    for i in range(1, m):
        thresholds[i] = first_u_reaching(patterns[i])
    thresholds[m] = two_k

    probs = {}
    for i in range(m):
        count = thresholds[i + 1] - thresholds[i]
        assert count >= 0, (i, patterns[i], thresholds[i], thresholds[i + 1])
        probs[patterns[i]] = Fraction(count, two_k)
    assert sum(probs.values()) == 1
    return probs, entries

# ---- binning + chi-square ----

def quantile_bins(entries, probs, num_bins):
    """Partition the value-sorted pattern list into ~num_bins contiguous
    groups of roughly equal expected mass (never splitting a pattern across
    a bin boundary). Returns a list of (patterns_in_bin, total_prob)."""
    bins = []
    cur_patterns = []
    cur_prob = Fraction(0)
    target = Fraction(1, num_bins)
    for _, pattern in entries:
        cur_patterns.append(pattern)
        cur_prob += probs[pattern]
        if cur_prob >= target and len(bins) < num_bins - 1:
            bins.append((cur_patterns, cur_prob))
            cur_patterns = []
            cur_prob = Fraction(0)
    if cur_patterns:
        bins.append((cur_patterns, cur_prob))
    return bins

def chi_square_stat(counts_by_pattern, bins, n_draws):
    stat = Fraction(0)
    rows = []
    for pats, prob in bins:
        observed = sum(counts_by_pattern.get(p, 0) for p in pats)
        expected = prob * n_draws
        term = (Fraction(observed) - expected) ** 2 / expected
        stat += term
        rows.append((observed, expected, term))
    return stat, rows

def parse_counts(path):
    """One raw pattern per line (hex 0x.. or decimal); returns pattern->count."""
    counts = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            p = int(line, 0)
            counts[p] = counts.get(p, 0) + 1
    return counts

# ---- CLI ----

def cmd_table(args):
    probs, entries = expected_probabilities(args.width)
    n = WIDTHS[args.width]
    print(f"{args.width}: {len(entries)} reachable patterns, k={KBITS[args.width]}")
    for v, p in entries:
        print(f"  0x{p:0{(n+3)//4}x}  value={float(v):.10g}  prob={probs[p]} ({float(probs[p]):.3e})")

def cmd_chi2(args):
    probs, entries = expected_probabilities(args.width)
    bins = quantile_bins(entries, probs, args.bins)
    counts = parse_counts(args.counts_file)
    n_draws = sum(counts.values())
    stat, rows = chi_square_stat(counts, bins, n_draws)
    dof = len(bins) - 1
    stat_f = float(stat)
    print(f"{args.width}: N={n_draws} draws, {len(bins)} bins, dof={dof}")
    print(f"chi-square statistic = {stat_f:.4f}")
    if _scipy_chi2 is not None:
        p_value = 1.0 - _scipy_chi2.cdf(stat_f, dof)
        print(f"p-value = {p_value:.4f}  (fail to reject H0 [uniform-on-[0,1)] at alpha=0.01 if p > 0.01)")
    else:
        print("(scipy not installed -- compare the statistic above to a "
              f"chi-square critical value table at dof={dof} yourself)")
    if args.verbose:
        for i, ((pats, prob), (observed, expected, term)) in enumerate(zip(bins, rows)):
            print(f"  bin {i:3d}: patterns={len(pats):5d} observed={observed:6d} "
                  f"expected={float(expected):9.2f} term={float(term):.4f}")

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)

    t = sub.add_parser('table', help='print the exact expected-probability table')
    t.add_argument('width', choices=WIDTHS.keys())
    t.set_defaults(func=cmd_table)

    c = sub.add_parser('chi2', help='chi-square ship-drawn counts against the exact table')
    c.add_argument('width', choices=['posit8', 'posit16'])
    c.add_argument('counts_file')
    c.add_argument('--bins', type=int, default=32)
    c.add_argument('--verbose', action='store_true')
    c.set_defaults(func=cmd_chi2)

    args = ap.parse_args()
    args.func(args)

if __name__ == '__main__':
    main()
