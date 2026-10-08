"""
Figure 1 panels for the calibrated model.
  A  cell amounts over 0-120 h (mean ± SD, ten held-out seeds)
  B  cytokine amounts over 0-120 h, each scaled to its maximum (mean ± SD)
  C  spatial snapshots (one simulation): debris field on the mesh with parcels drawn as dots coloured by
     their composition (blend of PMN blue, M0 grey, M1 red, M2 green, MSC yellow), the snapshot technique of
     the original Figure 1.
Usage: python tools/fig1_panels.py [--seed 11000] [--hours 0,24,48,72,96,120]
Writes out/figures/fig1/ (separate panels as PNG 600 dpi and PDF, plus an assembled preview).
"""
import argparse
import glob
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.lines import Line2D

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from utils import parse_nodes_elements            # noqa: E402
from domain_model import DomainModel              # noqa: E402
from element_agent_optimized import ElementAgent  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--seed", type=int, default=11000)
ap.add_argument("--hours", default="0,24,48,72,96,120")
ap.add_argument("--pmn_weight_late", type=float, default=0.1,
                help="weight of PMN in the colour blend of snapshots after 48 h (display only)")
args = ap.parse_args()
HOURS = [int(h) for h in args.hours.split(",")]
OUT = REPO / "out" / "figures" / "fig1"
OUT.mkdir(parents=True, exist_ok=True)

# colours of the original snapshot legend
BASE = {"PMN": (0, 0, 1), "M0": (0.5, 0.5, 0.5), "M1": (1, 0, 0), "M2": (0, 0.8, 0), "MSC": (1, 1, 0)}
LINE = {"PMN": "#0000ff", "M0": "#808080", "M1": "#ff0000", "M2": "#00cc00", "MSC": "#ffff00"}
CYTO = [("amt_c1", "TNF-α", "#7b3294"), ("amt_c2", "IL-10", "#008b8b"), ("amt_c3", "TGF-β", "#8c510a"), ("amt_c4", "VEGF-A", "#e66101")]
INK = "#1f1f1e"
plt.rcParams.update({"font.size": 8, "axes.edgecolor": INK, "axes.labelcolor": INK, "xtick.color": INK,
                     "ytick.color": INK, "axes.spines.top": False, "axes.spines.right": False})
YEL = [pe.withStroke(linewidth=2.4, foreground="#8a8a00")]   # outline so yellow is visible on white


def save(fig, name):
    fig.savefig(OUT / f"{name}.png", dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


# ------------------------------------------------------------------ snapshots (run one simulation)
snap_file = OUT / f"snapshots_seed{args.seed}.json"
nodes, elements = parse_nodes_elements(REPO / "data" / "node_elements.txt")
if not snap_file.exists():
    p = json.load(open(REPO / "params" / "model_of_record.json"))
    model = DomainModel(nodes, elements, p, seed=args.seed)
    model.enable_EC = False
    snaps = {}

    def grab(t):
        snaps[t] = {"debris": {str(e): float(v) for e, v in model.debris_field.items()},
                    "parcels": [[int(a.element_id) if str(a.element_id).isdigit() else a.element_id,
                                 *[float(a.state[i]) for i in (0, 1, 2, 3, 7)]]
                                for a in model.schedule.agents if isinstance(a, ElementAgent)]}
    if 0 in HOURS:
        grab(0)
    while model.time < max(HOURS):
        model.step()
        if model.time in HOURS:
            grab(model.time)
    json.dump(snaps, open(snap_file, "w"))
snaps = {int(k): v for k, v in json.load(open(snap_file)).items()}

centroid = {str(e): np.mean([nodes[n] for n in v], axis=0) for e, v in elements.items()}
polys = [[nodes[n] for n in elements[e][:4]] for e in elements]
keys = [str(e) for e in elements]
rng = np.random.default_rng(0)
cmap = plt.cm.viridis
norm = plt.Normalize(0, 200)


def snapshot(ax, t, cbar=False):
    s = snaps[t]
    vals = np.array([s["debris"][k] for k in keys])
    pc = PolyCollection(polys, array=vals, cmap=cmap, norm=norm, edgecolors="grey", linewidths=0.15, alpha=0.4)
    ax.add_collection(pc)
    xy, col, alp = [], [], []
    wN = args.pmn_weight_late if t > 48 else 1.0      # display only: PMN dimmed after 48 h
    for e, N, M0, M1, M2, C in s["parcels"]:
        amt = np.array([wN * N, M0, M1, M2, C])
        if amt.sum() <= 0:
            continue
        w = amt / amt.sum()
        col.append(np.clip(sum(wi * np.array(BASE[k]) for wi, k in zip(w, BASE)), 0, 1))
        alp.append(0.25 if (t > 48 and w[0] > 0.5) else 1.0)     # PMN-dominated parcels faded
        xy.append(centroid[str(e)] + rng.uniform(-0.05, 0.05, 2))
    if xy:
        xy, col, alp = np.array(xy), np.array(col), np.array(alp)
        order = np.argsort(alp)                             # faded PMN parcels drawn underneath
        ax.scatter(xy[order, 0], xy[order, 1], c=np.column_stack([col, alp])[order], s=6, linewidths=0.15,
                   edgecolors=[(0.2, 0.2, 0.2, a) for a in alp[order]], zorder=3)
    ax.set_xlim(-2.2, 7.7); ax.set_ylim(-10.3, 10.3); ax.set_aspect("equal")
    ax.set_title(f"{t} h", fontsize=9, color=INK)
    ax.set_xlabel("x (mm)")
    return pc


legend_handles = [Line2D([0], [0], marker="o", ls="none", markersize=6, markerfacecolor=BASE[k],
                         markeredgecolor="#333333", markeredgewidth=0.4, label=k) for k in BASE]

# single panels
for t in HOURS:
    fig, ax = plt.subplots(figsize=(1.9, 3.6))
    snapshot(ax, t)
    ax.set_ylabel("y (mm)")
    save(fig, f"C_snapshot_{t:03d}h")

# row of snapshots with legend and debris colour bar
fig, axes = plt.subplots(1, len(HOURS), figsize=(1.55 * len(HOURS) + 1.4, 3.6), sharey=True)
for i, (ax, t) in enumerate(zip(axes, HOURS)):
    pc = snapshot(ax, t)
    if i:
        ax.tick_params(labelleft=False)
axes[0].set_ylabel("y (mm)")
cb = fig.colorbar(pc, ax=axes, fraction=0.025, pad=0.02)
cb.set_label("debris (a.u.·mm⁻²)"); cb.outline.set_visible(False); cb.solids.set_alpha(0.4)
fig.legend(handles=legend_handles, loc="lower center", ncol=5, frameon=False, bbox_to_anchor=(0.45, -0.04),
           title="parcel colour = blend of its composition", title_fontsize=7.5, fontsize=7.5)
save(fig, "C_snapshots_row")

# ------------------------------------------------------------------ A and B (ten held-out seeds)
line_handles = [Line2D([0], [0], color=LINE[k], lw=1.8, label=k, path_effects=YEL if k == "MSC" else None) for k in BASE]
B = [json.load(open(f)) for f in sorted(glob.glob(str(REPO / "out" / "final" / "timecourse" / "baseline_*.json")))]
t = np.array([h["t"] for h in B[0]["hourly"]])


def series(k):
    Y = np.array([[h[k] for h in r["hourly"]] for r in B])
    return Y.mean(0), Y.std(0, ddof=1)


fig, ax = plt.subplots(figsize=(3.6, 2.6))
for k in BASE:
    m, s = series(k)
    ax.fill_between(t, m - s, m + s, color=LINE[k], alpha=0.2, lw=0)
    ln, = ax.plot(t, m, color=LINE[k], lw=1.6, label=k)
    if k == "MSC":
        ln.set_path_effects(YEL)
ax.set_xlim(0, 120); ax.set_xticks([0, 24, 48, 72, 96, 120])
ax.set_xlabel("time after osteotomy (h)"); ax.set_ylabel("cell amount (cell-equivalents)")
ax.legend(handles=line_handles, frameon=False, fontsize=7, ncol=1, loc="upper left")
ax.grid(axis="y", color="#e2e1dd", lw=0.6)
save(fig, "A_cells")

fig, ax = plt.subplots(figsize=(3.0, 2.6))
for k, lab, c in CYTO:
    m, s = series(k)
    mx = m.max()
    ax.fill_between(t, (m - s) / mx, (m + s) / mx, color=c, alpha=0.2, lw=0)
    ax.plot(t, m / mx, color=c, lw=1.6, label=lab)
ax.set_xlim(0, 120); ax.set_ylim(0, 1.12); ax.set_xticks([0, 24, 48, 72, 96, 120])
ax.set_xlabel("time after osteotomy (h)"); ax.set_ylabel("amount / maximum")
ax.legend(frameon=False, fontsize=7, loc="upper left")
ax.grid(axis="y", color="#e2e1dd", lw=0.6)
save(fig, "B_cytokines")

# ------------------------------------------------------------------ assembled preview (A, B over C)
fig = plt.figure(figsize=(7.2, 6.4))
outer = fig.add_gridspec(2, 1, height_ratios=[1, 1.55], hspace=0.28)
top = outer[0].subgridspec(1, 2, width_ratios=[1.35, 1], wspace=0.32)
bot = outer[1].subgridspec(1, len(HOURS), wspace=0.06)
axA, axB = fig.add_subplot(top[0]), fig.add_subplot(top[1])
for k in BASE:
    m, s_ = series(k)
    axA.fill_between(t, m - s_, m + s_, color=LINE[k], alpha=0.2, lw=0)
    ln, = axA.plot(t, m, color=LINE[k], lw=1.4)
    if k == "MSC":
        ln.set_path_effects(YEL)
for k, lab, c in CYTO:
    m, s_ = series(k)
    axB.fill_between(t, (m - s_) / m.max(), (m + s_) / m.max(), color=c, alpha=0.2, lw=0)
    axB.plot(t, m / m.max(), color=c, lw=1.4, label=lab)
axA.legend(handles=line_handles, frameon=False, fontsize=6.5, loc="upper left")
axB.legend(frameon=False, fontsize=6.5, loc="upper left")
for ax, yl in ((axA, "cell amount (cell-eq.)"), (axB, "amount / maximum")):
    ax.set_xlim(0, 120); ax.set_xticks([0, 24, 48, 72, 96, 120]); ax.set_xlabel("time (h)"); ax.set_ylabel(yl)
    ax.grid(axis="y", color="#e2e1dd", lw=0.6)
axB.set_ylim(0, 1.12)
axA.set_title("A", loc="left", fontweight="bold"); axB.set_title("B", loc="left", fontweight="bold")
axes = [fig.add_subplot(bot[i]) for i in range(len(HOURS))]
for i, (ax, h) in enumerate(zip(axes, HOURS)):
    pc = snapshot(ax, h)
    if i:
        ax.tick_params(labelleft=False)
    else:
        ax.set_ylabel("y (mm)")
axes[0].text(-0.45, 1.1, "C", transform=axes[0].transAxes, fontweight="bold", fontsize=9)
cb = fig.colorbar(pc, ax=axes, fraction=0.02, pad=0.015)
cb.set_label("debris (a.u.·mm⁻²)", fontsize=7); cb.outline.set_visible(False); cb.solids.set_alpha(0.4)
fig.legend(handles=legend_handles, loc="lower center", ncol=5, frameon=False, bbox_to_anchor=(0.47, 0.0),
           fontsize=7, title=f"parcel colour = blend of its composition (after 48 h, PMN weighted ×{args.pmn_weight_late:g}; "
                             "PMN-dominated parcels faded)", title_fontsize=7)
save(fig, "fig1_preview")
print("written to", OUT)
