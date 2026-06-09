#!/usr/bin/env python3
"""
make_pbe_figure.py — render a preliminary-data figure from pbe_pipeline output.

Usage:
    python3 make_pbe_figure.py structure.json  out.png  ["Title"]  [--representative]

Reads the structure.json written by pbe_pipeline.py and draws:
  - RBP-coverage track (partner-RBP footprint along the transcript)
  - accessibility track (per-base unpaired probability)
  - nominated PBE window(s) shaded across both tracks
  - candidate disruptor target positions as ticks
  - a small table of the top disruptors

Pass --representative to stamp an "ILLUSTRATIVE — replace with real data" banner.

NOTE: the figure is a labels/data panel intended as a separate supplemental attachment.
"""
import json, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

def main():
    if len(sys.argv) < 3:
        sys.exit("usage: make_pbe_figure.py structure.json out.png [Title] [--representative]")
    jpath, out = sys.argv[1], sys.argv[2]
    title = sys.argv[3] if len(sys.argv) > 3 and not sys.argv[3].startswith("--") else None
    rep = "--representative" in sys.argv

    d = json.load(open(jpath))
    n = d["length"]
    x = np.arange(n)
    cov = np.array(d["rbp_coverage"])
    unp = np.array(d["unpaired_prob"])
    pbe = d["pbe_windows"]
    disr = d.get("top_disruptors", [])[:8]
    name = d.get("name", "lncRNA")

    NAVY, BLUE, AMBER, TEAL, RED = "#1F3864", "#2E6DB4", "#E0922E", "#2E9B8F", "#C0504D"
    fig, axes = plt.subplots(2, 1, figsize=(11, 5.2), sharex=True,
                             gridspec_kw={"height_ratios": [1, 1], "hspace": 0.12})

    minw = max(1, int(n * 0.006))  # minimum visible width on a long transcript

    # PBE shading (min visible width)
    for w in pbe:
        ww = max(w["end"] - w["start"], minw)
        for ax in axes:
            ax.add_patch(Rectangle((w["start"], 0), ww, 1,
                                   transform=ax.get_xaxis_transform(),
                                   color=AMBER, alpha=0.30, zorder=0))

    # contiguous covered intervals (RBP footprint / motif sites), drawn with min width
    def intervals(arr):
        out, i = [], 0
        while i < len(arr):
            if arr[i] > 0:
                j = i
                while j < len(arr) and arr[j] > 0:
                    j += 1
                out.append((i, j)); i = j
            else:
                i += 1
        return out

    ax0 = axes[0]
    for (s, e) in intervals(cov):
        ax0.add_patch(Rectangle((s, 0), max(e - s, minw), 1.0,
                                transform=ax0.get_xaxis_transform(), color=RED, alpha=0.9, lw=0))
    ax0.set_ylabel("Partner-RBP\nfootprint", fontsize=10, color=NAVY)
    ax0.set_ylim(0, 1.05); ax0.set_yticks([0, 1])
    ax0.tick_params(labelsize=8)

    ax1 = axes[1]
    ax1.fill_between(x, unp, color=BLUE, alpha=0.75, lw=0)
    ax1.set_ylabel("Accessibility\n(unpaired prob.)", fontsize=10, color=NAVY)
    ax1.set_ylim(0, 1.05); ax1.set_yticks([0, 1])
    ax1.set_xlabel("Transcript position (nt)", fontsize=10, color=NAVY)
    ax1.tick_params(labelsize=8)

    # disruptor ticks on accessibility panel (min visible width)
    for c in disr:
        ax1.add_patch(Rectangle((c["target_start"], -0.16), max(c["len"], minw), 0.10,
                                transform=ax1.get_xaxis_transform(), clip_on=False,
                                color=TEAL, alpha=0.9))
    ax1.text(0, -0.22, "candidate disruptor target sites",
             transform=ax1.get_xaxis_transform(), fontsize=8, color=TEAL, va="top")

    # PBE label
    if pbe:
        w = pbe[0]
        ax0.annotate("PBE", (w["start"] + (w["end"] - w["start"]) / 2, 1.02),
                     xycoords=("data", "axes fraction"), ha="center", va="bottom",
                     fontsize=10, fontweight="bold", color="#8A5A12")

    ttl = title or f"In-silico PBE mapping and disruptor nomination — {name}"
    fig.suptitle(ttl, fontsize=12, fontweight="bold", color=NAVY, y=0.98)

    # top-disruptor table (compact, below)
    if disr:
        rows = [[f"{c['target_start']}-{c['target_end']}", c["disruptor_antisense_5to3"],
                 f"{c['Tm_C_approx']}", f"{c['GC_pct']}", f"{c['target_accessibility']}"]
                for c in disr[:4]]
        tbl = ax1.table(cellText=rows,
                        colLabels=["target", "disruptor 5'->3' (antisense)", "Tm~", "GC%", "access."],
                        cellLoc="center", colLoc="center",
                        bbox=[0.0, -0.72, 1.0, 0.40])
        tbl.auto_set_font_size(False); tbl.set_fontsize(7.5)
        for (r, _), cell in tbl.get_celld().items():
            cell.set_edgecolor("#CCCCCC")
            if r == 0:
                cell.set_facecolor("#DCE6F4"); cell.set_text_props(weight="bold", color=NAVY)

    if rep:
        fig.text(0.5, 0.5, "REPRESENTATIVE — replace with real data", fontsize=26,
                 color="red", alpha=0.18, ha="center", va="center", rotation=18, weight="bold")

    fig.subplots_adjust(bottom=0.42, top=0.9, left=0.11, right=0.97)
    fig.savefig(out, dpi=200, facecolor="white")
    print("wrote", out)

if __name__ == "__main__":
    main()
