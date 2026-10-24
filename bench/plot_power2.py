#!/usr/bin/env python3
"""
Plot fox vs fox_shared from bench/results_power2.csv.
  python3 bench/plot_power2.py [csv] [out.png]

Top row: median compute_ms per np, one panel per N (each panel has its own
y scale). Bottom row: fox_shared time / fox time; above 1 means fox_shared is
slower. np beyond the machine's hardware threads is shaded (oversubscribed).
"""
import csv
import os
import statistics
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CSV  = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "bench", "results_power2.csv")
OUT  = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "report", "power2_fox_vs_shared.png")
HW_THREADS = int(os.environ.get("HW_THREADS", os.cpu_count() or 16))

SURFACE   = "#fcfcfb"
INK       = "#0b0b0b"
INK_2     = "#52514e"
GRID      = "#e4e3df"
VARIANTS  = {  # fixed categorical order: slot 1 blue, slot 2 orange
    "fox":        dict(color="#2a78d6", marker="o", label="fox (bcast + shift)"),
    "fox_shared": dict(color="#eb6834", marker="s", label="fox_shared (local copy + allgather)"),
}

times = defaultdict(list)  # (n, variant, np) -> [ms]
with open(CSV) as f:
    for r in csv.DictReader(f):
        if r["compute_ms"] != "NA":
            times[(int(r["n"]), r["variant"], int(r["np"]))].append(float(r["compute_ms"]))

med = {k: statistics.median(v) for k, v in times.items()}
ns  = sorted({k[0] for k in med})
nps = sorted({k[2] for k in med})

plt.rcParams.update({
    "font.size": 9, "axes.edgecolor": INK_2, "axes.labelcolor": INK_2,
    "xtick.color": INK_2, "ytick.color": INK_2, "text.color": INK,
})
fig, axes = plt.subplots(2, len(ns), figsize=(3.0 * len(ns), 5.6),
                         sharex=True, facecolor=SURFACE,
                         gridspec_kw=dict(height_ratios=[3, 2]))

def style(ax):
    ax.set_facecolor(SURFACE)
    ax.set_xscale("log", base=2)
    ax.set_xticks(nps, [str(p) for p in nps])
    ax.minorticks_off()
    ax.grid(True, axis="y", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    over = [p for p in nps if p > HW_THREADS]
    if over:
        ax.axvspan(HW_THREADS * 1.4, max(nps) * 1.4, color=GRID, alpha=0.5, linewidth=0)

for col, n in enumerate(ns):
    top, bot = axes[0][col], axes[1][col]
    style(top)
    style(bot)

    for v, spec in VARIANTS.items():
        xs = [p for p in nps if (n, v, p) in med]
        top.plot(xs, [med[(n, v, p)] for p in xs], color=spec["color"],
                 marker=spec["marker"], markersize=6, linewidth=2,
                 markeredgecolor=SURFACE, markeredgewidth=1.5, label=spec["label"])
    top.set_title(f"N = {n}", color=INK, fontsize=10)
    top.set_ylim(bottom=0)

    xs = [p for p in nps if (n, "fox", p) in med and (n, "fox_shared", p) in med]
    ratio = [med[(n, "fox_shared", p)] / med[(n, "fox", p)] for p in xs]
    bot.axhline(1, color=INK_2, linewidth=1, linestyle="--")
    bot.plot(xs, ratio, color=VARIANTS["fox_shared"]["color"], marker="s",
             markersize=6, linewidth=2, markeredgecolor=SURFACE, markeredgewidth=1.5)
    for x, y in zip(xs, ratio):
        # below the point when under 1, so the label never sits on the reference line
        bot.annotate(f"{y:.2f}", (x, y), textcoords="offset points",
                     xytext=(0, 7 if y >= 1 else -12),
                     ha="center", fontsize=8, color=INK_2)
    bot.set_ylim(0, max(2.5, max(ratio) * 1.25))
    bot.set_xlabel("np")

axes[0][0].set_ylabel("compute_ms (median)")
axes[1][0].set_ylabel("fox_shared / fox")

handles, labels = axes[0][0].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper left", ncol=2, frameon=False,
           bbox_to_anchor=(0.01, 1.0))
fig.suptitle("Fox vs shared-copy Fox, min-plus APSP", x=0.99, ha="right",
             y=0.995, fontsize=11, color=INK)
fig.text(0.99, 0.005,
         f"shaded: np > {HW_THREADS} hardware threads (oversubscribed) · "
         "bottom row above 1 = fox_shared slower",
         ha="right", fontsize=8, color=INK_2)
fig.tight_layout(rect=(0, 0.02, 1, 0.95))
fig.savefig(OUT, dpi=200, facecolor=SURFACE)
print(OUT)
