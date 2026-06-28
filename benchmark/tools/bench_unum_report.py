#!/usr/bin/env python3
"""Assemble the /lib/unum benchmark table: interpreted Hoon vs Python/SoftUnum
vs jetted Hoon, per (posit width, arm).

Inputs: two scraped %bout grids (jetted at n=JN, interpreted at n=IN).  The
%bout line is `took <unit>/<dot-grouped>`; stripping the dots yields the time in
MICROSECONDS (the unit label is just magnitude).  Per-call = (arm - base) / n.
The Python column times SoftUnum (the C lib the jet calls) via ctypes, which is
the speed a Python user gets (call overhead included).

    /opt/anaconda3/bin/python tools/bench_unum_report.py JET.txt INTERP.txt JN IN
"""
import ctypes, os, re, sys, time

jet_f, int_f = sys.argv[1], sys.argv[2]
JN, IN = int(sys.argv[3]), int(sys.argv[4])

def parse(path):
    #  -> {(door,arm): microseconds}, plus ('fdp',door)
    out, label = {}, None
    for ln in open(path):
        m = re.search(r'\[%cell %(\w+) %([\w-]+)\]', ln)
        f = re.search(r'\[%fdp %(\w+)\]', ln)
        t = re.search(r'took\s+\S+?/([\d.]+)', ln)
        if m: label = (m.group(1), m.group(2))
        elif f: label = ('fdp', f.group(1))
        elif t and label:
            out[label] = int(t.group(1).replace('.', ''))   # dots -> µs
            label = None
    return out

jet, itp = parse(jet_f), parse(int_f)
DOORS = ['rpb', 'rph', 'rps']
ARMS = ['add','sub','mul','div','fma','sqt','neg','lth',
        'exp','log','sin','cos','atan','pow']

#  --- Python/SoftUnum (ctypes) per-call ---
HERE = os.path.dirname(os.path.abspath(__file__))
LIB = os.path.join(HERE, '..', '..', 'src', 'SoftPosit')  # placeholder
LIB = os.path.expanduser('~/urbit/SoftUnum/libsoftunum.dylib')
L = ctypes.CDLL(LIB)
u32 = ctypes.c_uint32
PFX = {'rpb': 'p8', 'rph': 'p16', 'rps': 'p32'}
ONE = {'rpb': 0x40, 'rph': 0x4000, 'rps': 0x40000000}
HALF = {'rpb': 0x38, 'rph': 0x3800, 'rps': 0x38000000}
UN = {'sqt','neg','exp','log','sin','cos','atan'}
BIN = {'add','sub','mul','div','pow','lth'}

def bind(pfx, arm):
    nm = {'sqt':'sqrt'}.get(arm, arm)
    f = getattr(L, f'{pfx}_{nm}')
    n = 1 if arm in UN else (3 if arm=='fma' else 2)
    f.argtypes = [u32]*n; f.restype = ctypes.c_int if arm=='lth' else u32
    return f

def py_us(door, arm):
    pfx = PFX[door]; f = bind(pfx, arm); x = HALF[door]; y = ONE[door]
    args = (x,) if arm in UN else ((x,y,y) if arm=='fma' else (x,y))
    M = 200_000
    for _ in range(2000): f(*args)
    t0 = time.perf_counter()
    for _ in range(M): f(*args)
    return (time.perf_counter()-t0)/M*1e6

def percall(d, door, arm, n):
    base = d.get((door,'base')); v = d.get((door,arm))
    if base is None or v is None: return None
    return max(v-base, 0)/n

print(f"\n/lib/unum jet benchmark  (jetted n={JN:,}; interpreted n={IN:,})")
print("per-call microseconds; speedup = interpreted / jetted\n")
hdr = f"{'width':5} {'arm':5} {'interp µs':>11} {'jetted µs':>10} {'py/SU µs':>9} {'jet vs interp':>13}"
print(hdr); print('-'*len(hdr))
for door in DOORS:
    for arm in ARMS:
        ji = percall(jet, door, arm, JN); ii = percall(itp, door, arm, IN)
        try: pu = py_us(door, arm)
        except Exception: pu = None
        if ji is None or ii is None: continue
        sp = (ii/ji) if ji else float('inf')
        print(f"{door:5} {arm:5} {ii:11.2f} {ji:10.3f} "
              f"{(f'{pu:.3f}' if pu is not None else '-'):>9} {sp:12.0f}x")
    #  fdp: per-element = took / n (no base)
    fj = jet.get(('fdp',door)); fi = itp.get(('fdp',door))
    if fj and fi:
        fjp, fip = fj/JN, fi/IN
        print(f"{door:5} {'fdp':5} {fip:11.2f} {fjp:10.3f} {'-':>9} {fip/fjp:12.0f}x")
    print()
