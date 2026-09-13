#!/usr/bin/env python3
"""Classify jet-vs-Hoon differences from two lagoon-kick --raw dumps of
hoon/float-diff.hoon. Exit 0 if every difference is the known overflow
divergence (jet ±MAX where the Hoon gives ±inf, urbit/urbit#7426), else 1."""
import re, sys
MAX = {'rh': (0x7bff, 0x7c00), 'rs': (0x7f7fffff, 0x7f800000), 'rd': (0x7fefffffffffffff, 0x7ff0000000000000),
       'rq': ((0x7ffe << 112) | ((1 << 112) - 1), 0x7fff << 112)}
SIGN = {'rh': 1 << 15, 'rs': 1 << 31, 'rd': 1 << 63, 'rq': 1 << 127}
def rows(s):
    toks = re.findall(r'\[|\]|0x[0-9a-f]+', s); st = []; root = None
    for t in toks:
        if t == '[': st.append([])
        elif t == ']':
            n = st.pop(); n = (n[0], n[1])
            if st: st[-1].append(n)
            else: root = n
        else:
            if st: st[-1].append(t)
            else: root = t
    out = []; n = root
    while isinstance(n, tuple):
        row, n = n; flat = []
        while isinstance(row, tuple): flat.append(row[0]); row = row[1]
        flat.append(row); out.append(flat)
    return out
def tas(x): return bytes.fromhex(x[2:].rjust(2, '0'))[::-1].decode()
J = rows(open(sys.argv[1]).read().strip().splitlines()[-1]); H = rows(open(sys.argv[2]).read().strip().splitlines()[-1])
assert len(J) == len(H), (len(J), len(H))
ops = ['add', 'sub', 'mul', 'div']; total = 0; overflow = 0; other = []
for j, h in zip(J, H):
    assert j[:4] == h[:4]
    d = tas(j[0]); m = tas(j[1])
    for k in range(4, 8):
        if j[k] == h[k]: continue
        total += 1
        jv, hv = int(j[k], 16), int(h[k], 16)
        mx, inf = MAX[d]
        if (jv & ~SIGN[d]) == mx and (hv & ~SIGN[d]) == inf and (jv & SIGN[d]) == (hv & SIGN[d]) and m != 'n':
            overflow += 1
        else:
            other.append(f"{d} %{m} {ops[k-4]} a={j[2]} b={j[3]} jet={j[k]} hoon={h[k]}")
print(f"{len(J)} rows, {total} differences: {overflow} overflow-to-inf (known), {len(other)} other")
for line in other[:40]: print("  ", line)
sys.exit(1 if other else 0)
