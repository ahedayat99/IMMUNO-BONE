"""
Supplementary figures from the calibrated-model outputs (ten held-out seeds).
  figS1_composition.png            parcel phenotype composition (M1 and M2 share) at 24, 72 and 120 h
  figS3_perturbation_cells.png     whole-callus cell amounts over 0-120 h, every perturbation scenario
  figS4_perturbation_cytokines.png domain cytokine amounts over 0-120 h, every perturbation scenario
Usage: python tools/supp_figs.py
"""
import glob
import json
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter

REPO = Path(__file__).resolve().parents[1]
FIN = REPO / "out" / "final"
OUT = REPO / "out" / "figures"
INK, GRID = "#1f1f1e", "#e4e3df"
plt.rcParams.update({"font.size": 8.5, "font.weight": "bold", "axes.labelweight": "bold", "axes.titleweight": "bold",
                     "axes.labelsize": 9, "xtick.labelsize": 8, "ytick.labelsize": 8, "axes.linewidth": 1.0,
                     "axes.edgecolor": INK, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
                     "grid.linewidth": 0.5, "font.family": "DejaVu Sans"})
POP = {"PMN": "#0000ff", "M0": "#808080", "M1": "#ff0000", "M2": "#00cc00", "MSC": "#ffff00"}
YEL = [pe.withStroke(linewidth=2.6, foreground="#8a8a00")]


def save(fig, name):
    fig.savefig(OUT / f"{name}.png", dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


runs = {}
for f in sorted(glob.glob(str(FIN / "timecourse" / "*.json"))):
    m = re.match(r"(.+)_(\d+)\.json", Path(f).name)
    runs.setdefault(m.group(1), []).append(json.load(open(f)))
t = np.arange(121)


def mean_series(name, key):
    return np.array([[h[key] for h in r["hourly"]] for r in runs[name]]).mean(0)


# ------------------------------------------------------------------ Fig. S1: parcel composition
fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.4), sharex=True)
bins = np.linspace(0, 1, 21)
for j, h in enumerate(("24", "72", "120")):
    comp = np.array([c for r in runs["baseline"] for c in r["composition"][h]])
    tot = comp.sum(1)
    for i, (k, col, idx) in enumerate((("M1", POP["M1"], 1), ("M2", POP["M2"], 2))):
        share = comp[:, idx] / tot
        w = np.full(len(share), 100 / len(share))
        ax = axes[i, j]
        ax.hist(share, bins=bins, weights=w, color=col, alpha=0.85, edgecolor="white", linewidth=0.8)
        p10, p50, p90 = np.percentile(share, (10, 50, 90))
        ax.axvline(p50, color=INK, lw=1.2, ls="--")
        ax.text(0.97, 0.95, f"median {100 * p50:.0f} %\n10–90 %: {100 * p10:.0f}–{100 * p90:.0f} %",
                transform=ax.transAxes, ha="right", va="top", fontsize=7.5, color=INK,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=1.5))
        if i == 0:
            ax.set_title(f"{h} h  (n = {len(share):,} parcels)", fontsize=9)
        if j == 0:
            ax.set_ylabel(f"parcels (%)\n{k} share")
        if i == 1:
            ax.set_xlabel(f"{k}/(M0+M1+M2) in parcel")
        else:
            ax.set_xlabel("M1/(M0+M1+M2) in parcel")
            ax.tick_params(labelbottom=True)
fig.tight_layout()
save(fig, "figS1_composition")

# ------------------------------------------------------------------ Figs. S3, S4: perturbation time courses
GROUPS = [
    ("Injury", [("debris_half", "debris ×0.5", "#1baf7a", "-"), ("debris_double", "debris ×2", "#eb6834", "-"),
                ("debris_quad", "debris ×4", "#eda100", "-")]),
    ("Progenitors", [("kpm_half", "k_pm ×0.5", "#2a78d6", "--"), ("kpm_double", "k_pm ×2", "#2a78d6", "-"),
                     ("klm_half", "K_lm ×0.5", "#eb6834", "--"), ("klm_double", "K_lm ×2", "#eb6834", "-"),
                     ("msc_half", "initial progenitors ×0.5", "#1baf7a", "--"), ("msc_double", "initial progenitors ×2", "#1baf7a", "-"),
                     ("kcm_half", "k_cm ×0.5", "#e87ba4", "--"), ("kcm_double", "k_cm ×2", "#e87ba4", "-")]),
    ("Polarization", [("m1_bias", "M1 bias", "#eb6834", "-"), ("m2_bias", "M2 bias", "#1baf7a", "-")]),
    ("Neutrophils", [("pmn_10pct", "initial PMN ×0.1", "#2a78d6", ":"), ("pmn_half", "initial PMN ×0.5", "#2a78d6", "--"),
                     ("pmn_double", "initial PMN ×2", "#2a78d6", "-"), ("krp_half", "k_rp ×0.5", "#eb6834", "--"),
                     ("krp_double", "k_rp ×2", "#eb6834", "-")]),
]


def grid(rows, fname, ylab):
    fig, axes = plt.subplots(len(rows), 4, figsize=(7.4, 1.55 * len(rows) + 1.3), sharex=True)
    for j, (gname, scen) in enumerate(GROUPS):
        for i, (key, lab) in enumerate(rows):
            ax = axes[i, j]
            ax.plot(t, mean_series("baseline", key), color=INK, lw=2.0, label="baseline")
            for name, slab, col, ls in scen:
                ax.plot(t, mean_series(name, key), color=col, lw=1.4, ls=ls, label=slab)
            ax.set_xticks([0, 48, 96]); ax.set_xlim(0, 120)
            ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}" if abs(v) >= 10 or v == 0 else f"{v:.2g}"))
            if j == 0:
                ax.set_ylabel(f"{lab}\n{ylab}")
            if i == 0:
                ax.set_title(gname, fontsize=9.5)
            if i == len(rows) - 1:
                ax.set_xlabel("time (h)")
        h, l = axes[0, j].get_legend_handles_labels()
        axes[-1, j].legend(h, l, frameon=False, fontsize=6.8, loc="upper center", bbox_to_anchor=(0.5, -0.38),
                           handlelength=2.2)
    fig.tight_layout(h_pad=0.6, w_pad=0.6)
    save(fig, fname)


grid([("PMN", "PMN"), ("M0", "M0"), ("M1", "M1"), ("M2", "M2"), ("MSC", "MSC")], "figS3_perturbation_cells", "(cell-eq.)")
grid([("amt_c1", "TNF-α"), ("amt_c2", "IL-10"), ("amt_c3", "TGF-β"), ("amt_c4", "VEGF-A")], "figS4_perturbation_cytokines", "(a.u.)")
print("written: figS1_composition, figS3_perturbation_cells, figS4_perturbation_cytokines")
