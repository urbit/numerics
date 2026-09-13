#!/usr/bin/env python3
"""Generate hoon/float-diff.hoon: for each door and mode, every (a, b) pair
over the safe+big value set with add/sub/mul/div results, as a flat list of
[mode a b add sub mul div]. Run it twice with lagoon-kick --raw, once with
the float jets (default) and once with HOON_FLOAT_JET_DISABLE=1, and diff:

  scripts/gen-float-diff.py
  (cd hoon && hoonc --arbitrary --output float-diff.jam float-diff.hoon .)
  ./target/release/lagoon-kick hoon/float-diff.jam --raw > /tmp/jet.txt
  HOON_FLOAT_JET_DISABLE=1 ./target/release/lagoon-kick hoon/float-diff.jam --raw > /tmp/hoon.txt
  scripts/float-diff-report.py /tmp/jet.txt /tmp/hoon.txt
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gen_float_values import doors
out = []
for d, (aura, safe, big) in doors.items():
    lit = ' '.join(f'`{aura}`{v}' for v in safe + big)
    out.append(f"=/  {d}s  `(list {aura})`~[{lit}]")
out.append(";:  weld")
for d, (aura, safe, big) in doors.items():
    for m in ('n', 'u', 'd', 'z'):
        out.append(f"  ^-  (list [@tas @tas @ @ @ @ @ @])")
        out.append(f"  %-  zing")
        out.append(f"  %+  turn  {d}s")
        out.append(f"  |=  a={aura}")
        out.append(f"  %+  turn  {d}s")
        out.append(f"  |=  b={aura}")
        out.append(f"  [%{d} %{m} a b (~(add {d} %{m}) a b) (~(sub {d} %{m}) a b) (~(mul {d} %{m}) a b) (~(div {d} %{m}) a b)]")
out.append("==")
path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'hoon', 'float-diff.hoon')
open(path, 'w').write('\n'.join(out) + '\n')
print('->', os.path.normpath(path))
