#!/usr/bin/env python3
"""
Plot seq vs fox from a run_bench.sh CSV.
  python3 bench/plot_lab.py [csv] [out.png] [title]

Left: minimum compute_ms per N, one line per np, seq as the reference.
Right: speedup against seq per N; above 1 means fox is faster.
"""
import csv
import os
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CSV  = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "bench", "report", "lab_rand0.csv")
OUT  = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "report", "lab_rand0_speedup.png")
TITLE = sys.argv[3] if len(sys.argv) > 3 else "Fox on LabCC, random matrices"

SURFACE = "#fcfcfb"
INK     = "#0b0b0b"
INK_2   = "#52514e"
GRID    = "#e4e3df"
SERIES  = [  # fixed categorical order: blue, orange, aqua, yellow, magenta
    dict(color="#2a78d6", marker="o"),
    dict(color="#eb6834", marker="s"),
    dict(color="#1baf7a", marker="^"),
    dict(color="#eda100", marker="P"),  # not "D": that is seq's marker
    dict(color="#e87ba4", marker="v"),
]

times = defaultdict(list)  # (n, impl, np) -> [ms]
with open(CSV) as f:
    for r in csv.DictReader(f):
        if r["compute_ms"] != "NA":
            times[(int(r["n"]), r["impl"], int(r["np"]))].append(float(r["compute_ms"]))

best = {k: min(v) for k, v in times.items()}
ns   = sorted({k[0] for k in best})
nps  = sorted({k[2] for k in best if k[1] == "fox"})
# np = 1 sits on the reference line by construction; its speedups stay unlabelled
labelled = [p for p in nps if p > 1]
seq  = {n: best[(n, "seq", 1)] for n in ns if (n, "seq", 1) in best}

plt.rcParams.update({
    "font.size": 9, "axes.edgecolor": INK_2, "axes.labelcolor": INK_2,
    "xtick.color": INK_2, "ytick.color": INK_2, "text.color": INK,
})
fig, (left, right) = plt.subplots(1, 2, figsize=(10, 4.2), facecolor=SURFACE)

def style(ax):
    ax.set_facecolor(SURFACE)
    ax.set_xscale("log", base=2)
    ax.set_xticks(ns, [str(n) for n in ns])
    ax.minorticks_off()
    ax.grid(True, axis="y", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.set_xlabel("N")

style(left)
style(right)

xs = sorted(seq)
# on top, so it stays visible where fox with np = 1 runs on the same line
left.plot(xs, [seq[n] for n in xs], color=INK_2, linestyle="--", linewidth=1.5,
          marker="D", markersize=5, markeredgecolor=SURFACE, label="seq", zorder=5)
if 1 in nps:
    # shown next to fox on one process: seq at speedup 1 for each N
    right.plot(xs, [1] * len(xs), color=INK_2, linestyle="--", linewidth=1.5,
               marker="D", markersize=5, markeredgecolor=SURFACE, zorder=5)
else:
    right.axhline(1, color=INK_2, linewidth=1, linestyle="--")
placed = defaultdict(list)  # N -> speedups already labelled there
for np_, spec in zip(nps, SERIES):
    xs = [n for n in ns if (n, "fox", np_) in best]
    line = dict(color=spec["color"], marker=spec["marker"], markersize=6, linewidth=2,
                markeredgecolor=SURFACE, markeredgewidth=1.5, label=f"fox, P = {np_}")
    left.plot(xs, [best[(n, "fox", np_)] for n in xs], **line)

    xs = [n for n in xs if n in seq]
    sp = [seq[n] / best[(n, "fox", np_)] for n in xs]
    right.plot(xs, sp, **line)
    for x, y in zip(xs, sp):
        if np_ not in labelled or y < 0.5:  # too close to the axis; the table has them
            continue
        # a label close to one already placed would overlap it; the table has it
        if any(abs(y - other) < 0.3 for other in placed[x]):
            continue
        # off the reference line and away from the other series' labels
        if np_ == labelled[0]:
            below = y < 1.2 or x == xs[-1]
            off = (0, -13 if below else 8)
        elif y < 1:  # would sit on the first series' labels; the table has it
            continue
        else:
            off = (-16, 6) if y < 1.5 or x == xs[-1] else (0, -13)
        right.annotate(f"{y:.2f}", (x, y), textcoords="offset points",
                       xytext=off, ha="center", fontsize=7.5, color=INK_2)
        placed[x].append(y)

left.set_yscale("log")
left.minorticks_off()
left.set_ylabel("compute_ms (min of runs, log)")
left.set_title("Execution time", color=INK, fontsize=10)

top = max(seq[n] / best[(n, "fox", p)] for (n, impl, p) in best if impl == "fox" and n in seq)
right.set_ylim(0, max(1.5, top * 1.2))
right.set_ylabel("speedup  T_seq / T_fox")
right.set_title("Speedup against seq", color=INK, fontsize=10)

handles, labels = left.get_legend_handles_labels()
rows = 1 if len(handles) <= 4 else 2  # more would run into the title
fig.legend(handles, labels, loc="upper left", ncol=-(-len(handles) // rows), frameon=False,
           bbox_to_anchor=(0.01, 1.0))
fig.suptitle(TITLE, x=0.99, ha="right",
             y=0.995, fontsize=11, color=INK)
fig.text(0.99, 0.005, "seq: speedup 1 on the right" if 1 in nps
         else "dashed line on the right: speedup 1 (fox as fast as seq)",
         ha="right", fontsize=8, color=INK_2)
fig.tight_layout(rect=(0, 0.02, 1, 0.95 if rows == 1 else 0.90))
fig.savefig(OUT, dpi=200, facecolor=SURFACE)
print(OUT)
