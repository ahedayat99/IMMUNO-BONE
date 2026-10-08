"""
Main data figures for the revised manuscript, from the model-of-record outputs (10 held-out seeds).
  fig1_baseline.png        A: cell populations, B: cytokine amounts (normalized), C: spatial M2 share maps
  fig3B_convergence.png    calibration convergence (differential evolution + refinement)
  fig3CD_fit.png           C: macrophage ratios, D: cytokine folds — model vs experiment
  fig4_perturbations.png   effect sizes at 120 h (% change vs baseline)
  fig5_debris.png          debris clearance across scenarios
  figX_regions.png         A: region map, B: periosteal/centre M2-fraction ratio over time vs experiment
Usage: python tools/make_figures.py
"""
import glob
import json
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.colors import LinearSegmentedColormap, SymLogNorm
from matplotlib.ticker import FuncFormatter

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from utils import parse_nodes_elements  # noqa: E402

OUT = REPO / "out" / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIN = REPO / "out" / "final"
INK, MUTED, GRID, SURF = "#1f1f1e", "#6b6b67", "#e4e3df", "#fcfcfb"
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]          # fixed categorical order
SEQ = LinearSegmentedColormap.from_list("seq", ["#c6dbf5", "#2a78d6", "#0d3b75"])   # single-hue sequential
DIV = LinearSegmentedColormap.from_list("div", ["#2a78d6", "#f0efec", "#e34948"])   # diverging pair
plt.rcParams.update({"font.size": 8.5, "font.weight": "bold", "axes.labelweight": "bold", "axes.titleweight": "bold",
                     "axes.labelsize": 9, "xtick.labelsize": 8, "ytick.labelsize": 8, "axes.linewidth": 1.0,
                     "axes.edgecolor": INK, "axes.labelcolor": INK, "xtick.color": INK,
                     "ytick.color": INK, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.5, "figure.dpi": 300,
                     "savefig.dpi": 300, "font.family": "DejaVu Sans"})


def load_runs():
    runs = {}
    for f in sorted(glob.glob(str(FIN / "timecourse" / "*.json"))):
        m = re.match(r"(.+)_(\d+)\.json", Path(f).name)
        runs.setdefault(m.group(1), []).append(json.load(open(f)))
    return runs


runs = load_runs()
B = runs["baseline"]
t = np.arange(121)


def series(R, key):
    a = np.array([[r["hourly"][i][key] for i in range(121)] for r in R])
    return a.mean(0), a.std(0, ddof=1)


def end_label(ax, x, y, text, color):
    ax.annotate(text, (x, y), xytext=(4, 0), textcoords="offset points", va="center", fontsize=8, color=INK)
    ax.plot([x], [y], "o", ms=3, color=color)


# ---------------------------------------------------------------- Figure 1
nodes, elements = parse_nodes_elements(REPO / "data" / "node_elements.txt")
polys = {str(e): [nodes[n] for n in v[:4]] for e, v in elements.items()}
fig = plt.figure(figsize=(7.2, 6.4))
gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.15], hspace=0.38, wspace=0.35)
ax = fig.add_subplot(gs[0, 0:2])
ends = []
for c, (k, lab) in zip(CAT, [("PMN", "PMN"), ("M0", "M0"), ("M1", "M1"), ("M2", "M2"), ("MSC", "MSC")]):
    m, s = series(B, k)
    ax.fill_between(t, m - s, m + s, color=c, alpha=0.18, lw=0)
    ax.plot(t, m, color=c, lw=1.6, label=lab)
    ax.plot([120], [m[-1]], "o", ms=3, color=c)
    ends.append([m[-1], lab])
ends.sort()
y0, y1 = ax.get_ylim()
gap = 0.055 * (y1 - y0)                             # keep end labels ~1 text line apart
for i in range(1, len(ends)):
    ends[i][0] = max(ends[i][0], ends[i - 1][0] + gap)
for y, lab in ends:
    ax.annotate(lab, (120, y), xytext=(4, 0), textcoords="offset points", va="center", fontsize=8, color=INK)
ax.set_xlim(0, 132); ax.set_xticks([0, 24, 48, 72, 96, 120])
ax.set_xlabel("time after osteotomy (h)"); ax.set_ylabel("cell amount (cell-equivalents)")
ax.legend(frameon=False, ncol=5, loc="upper left", fontsize=8)
ax.set_title("A  Cell populations (mean ± SD, 10 simulations)", loc="left", fontsize=9, color=INK)
ax2 = fig.add_subplot(gs[0, 2])
for c, (k, lab) in zip(CAT[:4], [("amt_c1", "TNF-α"), ("amt_c2", "IL-10"), ("amt_c3", "TGF-β"), ("amt_c4", "VEGF-A")]):
    m, s = series(B, k)
    mx = m.max()
    ax2.fill_between(t, (m - s) / mx, (m + s) / mx, color=c, alpha=0.18, lw=0)
    ax2.plot(t, m / mx, color=c, lw=1.6, label=lab)
ax2.set_xticks([0, 48, 96]); ax2.set_xlim(0, 120); ax2.set_ylim(0, 1.15)
ax2.set_xlabel("time (h)"); ax2.set_ylabel("amount / maximum")
ax2.legend(frameon=False, fontsize=7.5, loc="upper left")
ax2.set_title("B  Cytokines", loc="left", fontsize=9, color=INK)
snap = next(r for r in B if r["seed"] == 11000)["snapshots"]
for j, hh in enumerate(["24", "72", "120"]):
    a = fig.add_subplot(gs[1, j])
    S = snap[hh]
    verts, vals, empty = [], [], []
    for e, v in S.items():
        M0, M1, M2 = v["cells"][1:4]
        if M0 + M1 + M2 > 1e-3:
            verts.append(polys[e]); vals.append(M2 / (M0 + M1 + M2))
        else:
            empty.append(polys[e])
    a.add_collection(PolyCollection(empty, facecolors="#dcdbd6", edgecolors="none"))
    pc = PolyCollection(verts, array=np.array(vals), cmap=SEQ, clim=(0, 0.8), edgecolors="none")
    a.add_collection(pc)
    a.set_xlim(-2.1, 7.6); a.set_ylim(-10.2, 10.2); a.set_aspect("equal"); a.grid(False)
    a.set_xticks([0, 4]); a.set_yticks([-10, 0, 10]); a.set_xlabel("x (mm)")
    if j == 0:
        a.set_ylabel("y (mm)")
    a.set_title(("C  " if j == 0 else "") + f"{int(hh)} h", loc="left", fontsize=9, color=INK)
cb = fig.colorbar(pc, ax=fig.axes[-3:], fraction=0.03, pad=0.02)
cb.set_label("M2 share of macrophages, M2/(M0+M1+M2)", fontsize=8); cb.outline.set_visible(False)
fig.text(0.01, 0.005, "C: one representative simulation (seed 11000); grey = no macrophages present. Plate side at x < −2 mm (not meshed); "
         "the gap spans |y| ≤ 1 mm.", fontsize=7, color=MUTED)
fig.savefig(OUT / "fig1_baseline.png", bbox_inches="tight", facecolor="white")
plt.close(fig)

# ---------------------------------------------------------------- Figure 3B
cal = REPO / "out" / "final" / "calibration" / "out" / "joint" / "A2_P2"
h1 = json.load(open(cal / "de_checkpoint.json"))["history"]
h2 = json.load(open(cal / "polish" / "de_checkpoint.json"))["history"]
g = list(range(1, len(h1) + 1)) + list(range(len(h1) + 1, len(h1) + len(h2) + 1))
best = [h["best"] for h in h1] + [h["best"] for h in h2]
mean = [h["mean"] for h in h1] + [h["mean"] for h in h2]
fig, ax = plt.subplots(figsize=(3.5, 2.4))
ax.plot(g, mean, color=CAT[1], lw=1.2, label="population mean")
ax.plot(g, best, color=CAT[0], lw=1.8, label="best")
ax.axvline(len(h1) + 0.5, color=MUTED, lw=0.8, ls="--")
ax.text(len(h1) + 1, max(mean) * 2.2, "refinement\n(6 seeds)", fontsize=7.5, color=MUTED, va="top")
ax.text(2, max(mean) * 2.2, "differential\nevolution (3 seeds)", fontsize=7.5, color=MUTED, va="top")
ax.set_yscale("log"); ax.set_ylim(min(best) * 0.6, max(mean) * 3); ax.set_xlabel("generation"); ax.set_ylabel("calibration objective")
ax.legend(frameon=False, fontsize=7.5, loc="center right")
ax.set_title("B  Calibration convergence", loc="left", fontsize=9, color=INK)
fig.savefig(OUT / "fig3B_convergence.png", bbox_inches="tight", facecolor="white")
plt.close(fig)

# ---------------------------------------------------------------- Figure 3C/D
H = [json.loads(l) for l in open(FIN / "heldout_n1.jsonl")]
exp_m = {"M1_M0_72": (1.60, 0.66), "M2_M0_72": (0.14, 0.20), "M1_M0_120": (1.69, 0.67), "M2_M0_120": (1.37, 0.59)}
labs = ["M1/M0\nDay 3", "M2/M0\nDay 3", "M1/M0\nDay 5", "M2/M0\nDay 5"]
fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 2.6), gridspec_kw={"wspace": 0.35})
x = np.arange(4)
em = [exp_m[k][0] for k in exp_m]; es = [exp_m[k][1] for k in exp_m]
mm = [np.mean([r["calibration"][k] for r in H]) for k in exp_m]; ms = [np.std([r["calibration"][k] for r in H], ddof=1) for k in exp_m]
a.errorbar(x - 0.12, em, yerr=es, fmt="s", ms=5, color="#8a8984", capsize=2, lw=1, label="experiment (mean ± SD)")
a.errorbar(x + 0.12, mm, yerr=ms, fmt="o", ms=5, color=CAT[0], capsize=2, lw=1, label="model (mean ± SD, held-out seeds)")
for xi, v in zip(x, mm):
    a.text(xi + 0.22, v, f"{v:.2f}", fontsize=7.5, va="center", color=INK)
a.set_xticks(x); a.set_xticklabels(labs); a.set_ylabel("ratio"); a.set_ylim(0, 3.0)
a.legend(frameon=False, fontsize=7.5, loc="upper left")
a.set_title("C  Macrophage ratios", loc="left", fontsize=9, color=INK)
cyt = [("c1", "TNF-α", 12.85), ("c2", "IL-10", 6.20), ("c3", "TGF-β", 8.97), ("c4", "VEGF-A", 5.54)]
cm = [np.mean([r["cytokine_fold"][c] for r in H]) for c, _, _ in cyt]; cs = [np.std([r["cytokine_fold"][c] for r in H], ddof=1) for c, _, _ in cyt]
b.scatter(x - 0.12, [v for _, _, v in cyt], marker="s", s=25, color="#8a8984", label="experiment (dPCR)", zorder=3)
b.errorbar(x + 0.12, cm, yerr=cs, fmt="o", ms=5, color=CAT[0], capsize=2, lw=1, label="model (amount fold)")
b.axhline(1, color=MUTED, lw=0.8, ls="--")
for xi, v in zip(x, cm):
    b.text(xi + 0.22, v, f"{v:.2f}", fontsize=7.5, va="center", color=INK)
b.set_yscale("log"); b.set_ylim(0.5, 20); b.set_xticks(x); b.set_xticklabels([n for _, n, _ in cyt])
b.set_ylabel("Day 5 / Day 3 fold change")
b.legend(frameon=False, fontsize=7.5, loc="upper right")
b.set_title("D  Cytokine fold-changes", loc="left", fontsize=9, color=INK)
fig.text(0.01, -0.09, "IL-10, TGF-β and VEGF-A folds were secondary calibration targets; TNF-α entered only through a direction penalty.",
         fontsize=7, color=MUTED)
fig.savefig(OUT / "fig3CD_fit.png", bbox_inches="tight", facecolor="white")
plt.close(fig)

# ---------------------------------------------------------------- Figure 4
summ = json.load(open(FIN / "summary.json"))
eff = summ["effect_120h_percent"]
order = ["debris_half", "debris_double", "debris_quad", "kpm_half", "kpm_double", "klm_half", "klm_double",
         "msc_half", "msc_double", "kcm_half", "kcm_double", "m1_bias", "m2_bias",
         "pmn_10pct", "pmn_half", "pmn_double", "krp_half", "krp_double"]
groups = [("Injury", 3), ("Progenitors", 8), ("Polarization", 2), ("Neutrophils", 5)]
rowlab = {"debris_half": "debris ×0.5", "debris_double": "debris ×2", "debris_quad": "debris ×4", "kpm_half": "proliferation k_pm ×0.5",
          "kpm_double": "proliferation k_pm ×2", "klm_half": "capacity K_lm ×0.5", "klm_double": "capacity K_lm ×2",
          "msc_half": "initial progenitors ×0.5", "msc_double": "initial progenitors ×2", "kcm_half": "progenitor VEGF k_cm ×0.5",
          "kcm_double": "progenitor VEGF k_cm ×2", "m1_bias": "M1 bias (k01 ×2, k02 ×0.5)", "m2_bias": "M2 bias (k01 ×0.5, k02 ×2)",
          "pmn_10pct": "initial PMN ×0.1", "pmn_half": "initial PMN ×0.5", "pmn_double": "initial PMN ×2",
          "krp_half": "PMN recruitment k_rp ×0.5", "krp_double": "PMN recruitment k_rp ×2"}
cols = ["PMN", "M0", "M1", "M2", "MSC", "amt_c1", "amt_c2", "amt_c3", "amt_c4", "debris"]
clab = ["PMN", "M0", "M1", "M2", "MSC", "TNF-α", "IL-10", "TGF-β", "VEGF-A", "debris"]
V = np.array([[eff[s][c] for c in cols] for s in order])
L = np.log2(1 + V / 100)                            # colour = log2 fold vs baseline: x0.5 and x2 get equal weight
fig, ax = plt.subplots(figsize=(7.2, 5.2))
im = ax.imshow(np.clip(L, -3, 3), cmap=DIV, vmin=-3, vmax=3, aspect="auto")
ax.set_xticks(np.arange(-.5, len(cols)), minor=True); ax.set_yticks(np.arange(-.5, len(order)), minor=True)
ax.grid(which="minor", color="white", linewidth=1.5); ax.grid(which="major", visible=False); ax.tick_params(which="minor", length=0)
for i in range(len(order)):
    for j in range(len(cols)):
        ax.text(j, i, f"{V[i, j]:+.0f}", ha="center", va="center", fontsize=7.5, color=INK)
ax.set_xticks(range(len(cols))); ax.set_xticklabels(clab); ax.xaxis.tick_top(); ax.tick_params(length=0)
ax.set_yticks(range(len(order))); ax.set_yticklabels([rowlab[s] for s in order])
for s in ax.spines.values():
    s.set_visible(False)
y0 = 0
for name, n in groups:
    ax.axhline(y0 - 0.5, color=INK, lw=0.8)
    ax.text(len(cols) - 0.3, y0 + n / 2 - 0.5, name, rotation=270, va="center", ha="left", fontsize=8, color=MUTED)
    y0 += n
cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.06, ticks=[-3, -2, -1, 0, 1, 2, 3])
cb.ax.set_yticklabels(["≤×1/8 (−88 %)", "×1/4 (−75 %)", "×1/2 (−50 %)", "×1 (0 %)", "×2 (+100 %)", "×4 (+300 %)", "≥×8 (+700 %)"],
                      fontsize=7.5); cb.outline.set_visible(False)
cb.set_label("fold change at 120 h vs baseline (log scale)", fontsize=8)
ax.text(0, -0.025, "Cell values: % change at 120 h vs baseline, mean of 10 simulations with the same random seeds as the baseline. Cytokines are domain amounts.",
         fontsize=7, color=MUTED, transform=ax.transAxes, va="top")
fig.savefig(OUT / "fig4_perturbations.png", bbox_inches="tight", facecolor="white")
plt.close(fig)

# ---------------------------------------------------------------- Figure 5
fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 2.6), gridspec_kw={"wspace": 0.3})
for ax_, log in ((a, True), (b, False)):
    for s in order:
        if s.startswith("debris"):
            continue
        m, _ = series(runs[s], "debris")
        ax_.plot(t, m, color="#c9c8c3", lw=0.7)
    for c, s, lab in ((CAT[2], "debris_half", "debris ×0.5"), (CAT[0], "baseline", "baseline"),
                      (CAT[1], "debris_double", "debris ×2"), (CAT[3], "debris_quad", "debris ×4")):
        m, _ = series(runs[s], "debris")
        ax_.plot(t, m, color=c, lw=1.6, label=lab)
    ax_.set_xticks([0, 24, 48, 72, 96, 120]); ax_.set_xlabel("time (h)")
    if log:
        ax_.set_yscale("log"); ax_.set_xlim(0, 124); ax_.set_ylabel("total debris (model units)")
    else:
        ax_.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v / 1e5:.0f}"))
        ax_.set_ylabel("total debris (×10⁵ model units)")
a.set_title("A  Debris clearance (log scale)", loc="left", fontsize=9, color=INK)
b.set_title("B  Linear scale", loc="left", fontsize=9, color=INK)
from matplotlib.lines import Line2D as _L2
h_, l_ = a.get_legend_handles_labels()
a.legend(h_ + [_L2([0], [0], color="#c9c8c3", lw=1)], l_ + ["other scenarios"], frameon=False, fontsize=7.5, loc="upper right")
fig.savefig(OUT / "fig5_debris.png", bbox_inches="tight", facecolor="white")
plt.close(fig)

# ---------------------------------------------------------------- New regional figure
fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.0), gridspec_kw={"width_ratios": [1, 2.2], "wspace": 0.3})
G = 1.88
zc, zp, zo = [], [], []
for e, v in elements.items():
    xy = np.mean([nodes[n] for n in v[:4]], axis=0)
    (zc if (xy[0] < G and abs(xy[1]) <= 1) else zp if (G <= xy[0] <= G + 2 and 1 <= abs(xy[1]) <= 8) else zo).append([nodes[n] for n in v[:4]])
a.add_collection(PolyCollection(zo, facecolors="#f0efec", edgecolors="none"))
a.add_collection(PolyCollection(zc, facecolors=CAT[1], edgecolors="none"))
a.add_collection(PolyCollection(zp, facecolors=CAT[0], edgecolors="none"))
a.set_xlim(-2.1, 7.6); a.set_ylim(-10.2, 10.2); a.set_aspect("equal"); a.grid(False)
a.set_xlabel("x (mm)"); a.set_ylabel("y (mm)")
a.text(-1.9, -1.4, "fracture\ncentre", color=INK, fontsize=7.5, va="top"); a.text(4.1, 6.5, "periosteal", color=INK, fontsize=7.5)
a.text(4.1, -7.0, "periosteal", color=INK, fontsize=7.5)
a.set_title("A  Regions", loc="left", fontsize=9, color=INK)
hours = sorted(int(k) for k in B[0]["zones"] if int(k) > 0)
def frac(z, hh, R):
    s = np.sum([r["zones"][str(hh)].get(z, [0] * 5) for r in R], axis=0)
    return s[3] / (s[1] + s[3])
pooled = [frac("periosteal", hh, B) / frac("fracture_centre", hh, B) for hh in hours]
for r in B:
    ps = []
    for hh in hours:
        zc_, zp_ = r["zones"][str(hh)].get("fracture_centre"), r["zones"][str(hh)].get("periosteal")
        ok = zc_ and zp_ and zc_[1] + zc_[3] > 0 and zp_[1] + zp_[3] > 0 and zc_[3] > 0
        ps.append((zp_[3] / (zp_[1] + zp_[3])) / (zc_[3] / (zc_[1] + zc_[3])) if ok else np.nan)
    b.plot(hours, ps, color="#c9c8c3", lw=0.7)
b.plot(hours, pooled, color=CAT[0], lw=1.8, marker="o", ms=3, label="model (pooled over 10 simulations; grey: individual)")
b.scatter([72, 120], [1.551, 1.948], marker="s", s=36, color=CAT[1], zorder=4, label="experiment (pooled cell counts, n = 4 rats)")
b.axhline(1, color=MUTED, lw=0.8, ls="--")
b.set_yscale("log"); b.set_ylim(0.2, 10); b.set_xticks([0, 24, 48, 72, 96, 120]); b.set_xlim(0, 124)
b.set_xlabel("time after osteotomy (h)"); b.set_ylabel("periosteal / centre M2 fraction, M2/(M0+M2)")
b.legend(frameon=True, facecolor="white", edgecolor="none", framealpha=0.9, fontsize=7.5, loc="upper right")
b.set_title("B  Periosteal M2 enrichment (not used in calibration)", loc="left", fontsize=9, color=INK)
fig.savefig(OUT / "figX_regions.png", bbox_inches="tight", facecolor="white")
plt.close(fig)
print("figures written to", OUT)
