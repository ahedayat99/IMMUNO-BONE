"""
Figure 3A (and Supplementary full version): one-at-a-time sensitivity analysis at the
pre-calibration (literature) base point. Signed normalized sensitivities, parameters ranked by
how many of the 8 outputs they rank top-10 in; the 14 selected parameters are marked.
Usage: python tools/fig_sa.py   (reads out/sa/sa_summary.json, writes out/figures/)
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams.update({"font.weight": "bold", "axes.labelweight": "bold", "axes.titleweight": "bold"})
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

REPO = Path(__file__).resolve().parents[1]
s = json.load(open(REPO / "out" / "sa" / "sa_summary.json"))
S, freq, selected = s["sensitivities"], s["frequency"], set(s["selected"])
OUTS = ["M1_M0_72", "M2_M0_72", "M1_M0_120", "M2_M0_120", "prod_c1", "prod_c2", "prod_c3", "prod_c4"]
COLS = ["M1/M0\nDay 3", "M2/M0\nDay 3", "M1/M0\nDay 5", "M2/M0\nDay 5", "TNF-α\nfold", "IL-10\nfold", "TGF-β\nfold", "VEGF-A\nfold"]
INK, MUTED, GRID = "#1f1f1e", "#6b6b67", "#d9d8d4"
CMAP = LinearSegmentedColormap.from_list("div", ["#2a78d6", "#f0efec", "#e34948"])  # diverging pair, grey midpoint
LIM = 0.6                                                                              # colour clipped at ±0.6
LABEL = {"k_e_pmn": "k_e,PMN", "a_M1_to_M2": "a_M1→M2", "a_M2_to_M1": "a_M2→M1", "d_c1": "d_c1", "d_c2": "d_c2",
         "d_c3": "d_c3", "d_c4": "d_c4", "d_m_p": "d_PMN", "d_c4_ec": "d_c4,EC", "k_e0": "k_e0", "k_e1": "k_e1", "k_e2": "k_e2"}


def order(names):
    return sorted(names, key=lambda k: (-freq[k], -max(abs(S[k][j]) for j in OUTS)))


def draw(names, path, title):
    V = np.array([[S[k][j] for j in OUTS] for k in names])
    h = 0.32 * len(names) + 1.6
    fig, ax = plt.subplots(figsize=(7.2, h), dpi=300)
    im = ax.imshow(np.clip(V, -LIM, LIM), cmap=CMAP, norm=TwoSlopeNorm(0, -LIM, LIM), aspect="auto")
    # 2-px surface gaps between cells
    ax.set_xticks(np.arange(-.5, len(OUTS)), minor=True)
    ax.set_yticks(np.arange(-.5, len(names)), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.5)
    ax.tick_params(which="minor", length=0)
    for i, k in enumerate(names):
        for j in range(len(OUTS)):
            v = V[i, j]
            if abs(v) >= 0.05:                         # selective labels: only non-negligible values
                ax.text(j, i, f"{v:+.2f}", ha="center", va="center", fontsize=7.5, color=INK)
    ax.set_xticks(range(len(OUTS)))
    ax.set_xticklabels(COLS, fontsize=8.5, color=INK)
    ax.xaxis.tick_top()
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels([("● " if k in selected else "   ") + LABEL.get(k, k) for k in names], fontsize=8.5, color=INK)
    for lab, k in zip(ax.get_yticklabels(), names):
        if k in selected:
            lab.set_fontweight("bold")
    ax.tick_params(length=0)
    for side in ax.spines.values():
        side.set_visible(False)
    # selection frequency column
    for i, k in enumerate(names):
        ax.text(len(OUTS) - 0.2, i, f"{freq[k]}/8", ha="left", va="center", fontsize=8.5,
                color=INK if k in selected else MUTED, fontweight="bold" if k in selected else "normal")
    ax.text(len(OUTS) - 0.2, -0.9, "top-10\nin", ha="left", va="bottom", fontsize=8, color=MUTED)
    cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.12, ticks=[-LIM, -0.3, 0, 0.3, LIM])
    cb.ax.set_yticklabels([f"≤−{LIM}", "−0.3", "0", "+0.3", f"≥+{LIM}"], fontsize=8, color=INK)
    cb.outline.set_visible(False)
    cb.set_label("normalized sensitivity S", fontsize=8.5, color=INK)
    ax.set_title(title, fontsize=8, color=INK, loc="left", pad=28)
    fig.text(0.01, 0.005, "● = selected for calibration (top 10 in ≥ 3 of 8 outputs). Cytokine outputs are Day 5/Day 3 production-rate folds. "
             "Values |S| < 0.05 not printed.", fontsize=7, color=MUTED)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


out = REPO / "out" / "figures"
out.mkdir(parents=True, exist_ok=True)
ranked = order(list(S))
draw(ranked[:20], out / "fig3A_sensitivity_top20.png", "One-at-a-time sensitivity at the literature base point (top 20 of 53 parameters)")
draw(ranked, out / "figS_sensitivity_all53.png", "One-at-a-time sensitivity at the literature base point (all 53 parameters)")
print("selected in top 20:", sum(k in selected for k in ranked[:20]), "of", len(selected))
print("ranked top 20:", ranked[:20])
