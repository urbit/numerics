#!/usr/bin/env python3
"""bench_plots.py — Tufte-style figures from bench-timing-all.csv.

Minimal ink: no chartjunk, range frames, log scale (data spans 7 orders of
magnitude), small multiples, direct labels.  Writes PNG+PDF to results/<date>/figures/.
"""
import csv, statistics, os, sys
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CSV = sys.argv[1] if len(sys.argv) > 1 else "results/2026-06-26/bench-timing-all.csv"
OUT = os.path.join(os.path.dirname(CSV), "figures"); os.makedirs(OUT, exist_ok=True)

# ---- load: per-call ns (base-subtracted), 32-bit ---------------------------
raw = defaultdict(dict)
for r in csv.DictReader(open(CSV)):
    if r["status"] != "ok" or r["word_size"] != "32":
        continue
    raw.setdefault((r["impl"], r["door"], r["arm"], int(r["n"])), {})[int(r["set"])] = float(r["elapsed_us"])

def percall(impl, door, arm):
    for n in (100000, 10000, 1000, 10):
        s = raw.get((impl, door, arm, n))
        if s:
            v = [us * 1000 / n for k, us in s.items() if k != 0]
            if v: return statistics.mean(v)
    return None

def pc_sub(impl, door, arm):
    a, b = percall(impl, door, arm), percall(impl, door, "base")
    return None if a is None else max(a - (b or 0), 1.0)

ARMS = ["exp","log","sin","cos","tan","atan","atan2","asin","acos","sqt","cbt","pow","pow-n","log-2","log-10"]
IMPLS = [("jetted", "Jetted (C)"), ("cheb", "Chebyshev, interpreted"), ("taylor", "Taylor, interpreted")]

# ---- Tufte rcParams --------------------------------------------------------
plt.rcParams.update({
    "font.family": "serif", "font.size": 9, "axes.titlesize": 10,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 3, "ytick.major.size": 0, "figure.dpi": 160,
    "axes.edgecolor": "#444444", "text.color": "#222222",
    "axes.labelcolor": "#222222", "xtick.color": "#444444", "ytick.color": "#222222",
})
INK = "#222222"; MID = "#888888"; LITE = "#bbbbbb"
MARK = {"jetted": dict(marker="o", ms=4.2, mfc=INK, mec=INK),
        "cheb":   dict(marker="o", ms=4.2, mfc="white", mec=INK, mew=0.9),
        "taylor": dict(marker="|", ms=8, mec=INK, mew=1.2)}

def rangeframe(ax, xs):
    lo, hi = min(xs), max(xs)
    ax.spines["bottom"].set_bounds(lo, hi)

def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(f"{OUT}/{name}.{ext}", bbox_inches="tight")
    plt.close(fig)

# === FIG 1: Cleveland dot plot — @rd, all arms, three implementations ========
order = sorted(ARMS, key=lambda a: (pc_sub("cheb","rd",a) or 0))
fig, ax = plt.subplots(figsize=(6.4, 4.4))
xs_all = []
for y, arm in enumerate(order):
    pts = [(impl, pc_sub(impl, "rd", arm)) for impl, _ in IMPLS]
    pts = [(i, v) for i, v in pts if v]
    if len(pts) > 1:  # connect with a thin guide line
        ax.plot([v for _, v in pts], [y]*len(pts), "-", color=LITE, lw=0.7, zorder=1)
    for impl, v in pts:
        ax.plot(v, y, zorder=3, **MARK[impl]); xs_all.append(v)
ax.set_yticks(range(len(order))); ax.set_yticklabels(order)
ax.set_xscale("log")
ax.set_xlabel("per-call time  (nanoseconds, log scale)")
ax.set_title("Math transcendentals on f64 (@rd): the jet vs interpreted Hoon",
             loc="left", pad=26)
rangeframe(ax, xs_all)
for x in (1e3,1e4,1e5,1e6,1e7,1e8):
    ax.axvline(x, color="#eeeeee", lw=0.6, zorder=0)
from matplotlib.lines import Line2D
handles = [Line2D([0],[0], ls="", **MARK[i], label=l) for i, l in IMPLS]
ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 1.09),
          ncol=3, frameon=False, fontsize=8, handletextpad=0.4, columnspacing=2.4)
save(fig, "fig1-rd-dotplot")

# === FIG 2: small multiples — jetted vs cheb-interp exp across all 4 doors ====
fig, axs = plt.subplots(1, 4, figsize=(7.6, 1.9), sharex=True)
for ax, door in zip(axs, ["rd","rs","rh","rq"]):
    j, c = pc_sub("jetted",door,"exp"), pc_sub("cheb",door,"exp")
    ax.plot([j, c], [0, 0], "-", color=LITE, lw=0.8)
    ax.plot(j, 0, **MARK["jetted"]); ax.plot(c, 0, **MARK["cheb"])
    ax.set_xscale("log"); ax.set_ylim(-1, 1); ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_title(f"@{door}   {c/j:.0f}×", fontsize=9, pad=6)
    ax.annotate(f"{j/1000:.1f}µs", (j,0), xytext=(0,-13), textcoords="offset points",
                ha="center", fontsize=7, color=INK)
    ax.annotate(f"{c/1000:.0f}µs", (c,0), xytext=(0,7), textcoords="offset points",
                ha="center", fontsize=7, color=INK)
axs[0].set_ylabel("exp", rotation=0, ha="right", va="center")
fig.suptitle("Jet speedup on exp is uniform across precisions   (jetted ●  →  Chebyshev ○, interpreted)",
             x=0.5, y=1.16, fontsize=10)
fig.text(0.5, -0.13, "per-call time (log scale)", ha="center", fontsize=8)
save(fig, "fig2-doors-exp")

# === FIG 3: accuracy — Taylor relative error vs Cheb (≤1 ULP), @rd ===========
acc = {}
try:
    for r in csv.DictReader(open(os.path.join(os.path.dirname(CSV), "bench-accuracy-rd.csv"))):
        acc[r["arm"]] = float(r["taylor_rel_err"])
except FileNotFoundError:
    acc = {}
if acc:
    eps = 2.0**-52
    order2 = sorted(acc, key=lambda a: acc[a])
    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    xs = []
    for y, arm in enumerate(order2):
        v = max(acc[arm], eps/2)
        ax.plot([eps, v], [y, y], "-", color=LITE, lw=0.7, zorder=1)  # from 1-ULP ref
        ax.plot(v, y, "o", ms=4.2, mfc=INK, mec=INK, zorder=3); xs.append(v)
    ax.axvline(eps, color=MID, lw=0.8, ls=(0,(4,3)))
    ax.annotate("1 ULP\n(Cheb/jet)", (eps, len(order2)-1), xytext=(4,0),
                textcoords="offset points", fontsize=7, color=MID, va="center")
    ax.set_yticks(range(len(order2))); ax.set_yticklabels(order2)
    ax.set_xscale("log"); ax.set_xlabel("Taylor relative error  (log scale)")
    ax.set_title("And the jet (Chebyshev) is more ACCURATE than legacy Taylor (@rd)", loc="left")
    rangeframe(ax, xs + [eps])
    save(fig, "fig3-accuracy")

print("wrote:", ", ".join(sorted(os.listdir(OUT))))
