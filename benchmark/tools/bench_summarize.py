#!/usr/bin/env python3
"""bench_summarize.py — aggregate bench-timing.csv into per-cell mean +/- stddev.

Input CSV columns (from bench_run.sh):
    word_size,impl,door,arm,n,set,elapsed_us,status

Per-call cost = elapsed_us * 1000 / n  (nanoseconds), with the matching %base
cell (same word_size,impl,door,n) subtracted to remove loop + fold overhead.
Set 0 is dropped as warm-up.  Failed cells (status != ok) are ignored.

Usage:
    tools/bench_summarize.py results/<date>/bench-timing.csv > results/<date>/bench-summary.tsv
"""
import csv
import sys
import statistics
from collections import defaultdict


def main(path):
    # rows[(ws,impl,door,arm,n)][set] = elapsed_us
    rows = defaultdict(dict)
    for r in csv.DictReader(open(path)):
        if r["status"] != "ok":
            continue
        key = (r["word_size"], r["impl"], r["door"], r["arm"], int(r["n"]))
        rows[key][int(r["set"])] = float(r["elapsed_us"])

    # per-call ns per cell, set 0 dropped
    def percall_ns(ws, impl, door, arm, n):
        sets = rows.get((ws, impl, door, arm, n), {})
        return {s: us * 1000.0 / n for s, us in sets.items() if s != 0}

    # baseline per (ws,impl,door,n): mean of %base per-call
    base_mean = {}
    for (ws, impl, door, arm, n) in rows:
        if arm == "base":
            pc = percall_ns(ws, impl, door, "base", n)
            if pc:
                base_mean[(ws, impl, door, n)] = statistics.mean(pc.values())

    out = csv.writer(sys.stdout, delimiter="\t")
    out.writerow(["word_size", "impl", "door", "arm", "n", "sets",
                  "mean_ns", "stddev_ns", "base_ns"])

    for (ws, impl, door, arm, n) in sorted(rows):
        if arm == "base":
            continue
        pc = percall_ns(ws, impl, door, arm, n)
        if not pc:
            continue
        base = base_mean.get((ws, impl, door, n), 0.0)
        vals = [v - base for v in pc.values()]
        mean = statistics.mean(vals)
        sd = statistics.pstdev(vals) if len(vals) > 1 else 0.0
        out.writerow([ws, impl, door, arm, n, len(vals),
                      f"{mean:.1f}", f"{sd:.1f}", f"{base:.1f}"])


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
