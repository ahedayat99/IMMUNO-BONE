"""
Figure 6: schematic of the IMMUNO-BONE coupling architecture (one 1-h step).
Usage: python tools/fig6_schematic.py   -> out/figures/fig6_schematic.png / .pdf
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Polygon, Wedge

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "out" / "figures"
INK, MUTED = "#1f1f1e", "#6b6b67"
POP = {"PMN": "#0000ff", "M0": "#808080", "M1": "#ff0000", "M2": "#00cc00", "MSC": "#ffff00"}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.weight": "bold"})

fig, ax = plt.subplots(figsize=(9.2, 5.4))
ax.set_xlim(0, 184); ax.set_ylim(0, 104); ax.axis("off")
BODY = 7.0


def box(x, y, w, h, fc, ec, title, body, num=None, ls="-"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.6,rounding_size=2.5", fc=fc, ec=ec, lw=1.6, ls=ls))
    tx = x + (7.0 if num else 2.0)
    if num:
        ax.add_patch(plt.Circle((x + 3.4, y + h - 4.0), 2.5, color=ec))
        ax.text(x + 3.4, y + h - 4.0, num, ha="center", va="center", color="white", fontsize=9, fontweight="bold")
    ax.text(tx, y + h - 4.0, title, ha="left", va="center", fontsize=9.5 if h > 20 else 8.5, color=INK, fontweight="bold")
    ax.text(x + 2.0, y + h - (8.8 if h > 20 else 7.6), body, ha="left", va="top", fontsize=BODY, color=INK, fontweight="normal", linespacing=1.4)


def arrow(p, q, rad=0.0, color=INK, ls="-", lw=1.8):
    ax.add_patch(FancyArrowPatch(p, q, connectionstyle=f"arc3,rad={rad}", arrowstyle="-|>", mutation_scale=14,
                                 color=color, lw=lw, ls=ls))


def pie(cx, cy, r, fracs):
    a = 90
    for k, f in fracs.items():
        if f > 0:
            ax.add_patch(Wedge((cx, cy), r, a, a + 360 * f, fc=POP[k], ec="#333333", lw=0.5))
            a += 360 * f


# ---------------------------------------------------------------- left: one element with mixed parcels
ax.text(2, 100, "One finite element e", fontsize=9.5, color=INK, fontweight="bold")
ax.add_patch(Polygon([(5, 57), (44, 60), (46, 91), (7, 89)], closed=True, fc="#eef3fb", ec="#2a78d6", lw=1.6))
pie(16, 79, 5.4, {"PMN": 0.15, "M0": 0.2, "M1": 0.45, "M2": 0.15, "MSC": 0.05})
pie(33, 81, 4.6, {"PMN": 0.6, "M0": 0.25, "M1": 0.15})
pie(26, 68, 5.0, {"M0": 0.15, "M1": 0.2, "M2": 0.45, "MSC": 0.2})
ax.text(26, 60.6, "≤ 3 parcels per element", ha="center", fontsize=6.8, color=INK, fontweight="normal")
ax.text(2, 51, "Parcel state", fontsize=8, color=INK)
ax.text(2, 47, "x_p = (N, M⁰, M¹, M², C)\ncontinuous amounts (cell-eq.)", fontsize=BODY, color=INK,
        fontweight="normal", va="top", linespacing=1.4)
ax.text(2, 36, "Element fields", fontsize=8, color=INK)
ax.text(2, 32, "debris D(e); cytokines c₁–c₄\n(TNF-α, IL-10, TGF-β, VEGF-A),\nshared by the resident parcels",
        fontsize=BODY, color=INK, fontweight="normal", va="top", linespacing=1.4)
for i, (k, c) in enumerate(POP.items()):
    ax.add_patch(plt.Circle((3 + i * 11.0, 12), 1.6, fc=c, ec="#333333", lw=0.5))
    ax.text(5.2 + i * 11.0, 12, k, va="center", fontsize=6.6, color=INK)
ax.text(2, 6, "phenotype = composition of the parcel", fontsize=6.8, color=MUTED, fontweight="normal")

# ---------------------------------------------------------------- right: the 1-h loop
C1, C2, C3, C4 = "#2a78d6", "#1baf7a", "#eb6834", "#b8860b"
box(56, 60, 58, 36, "#eaf2fc", C1, "Element reaction",
    "Coupled ODEs of all resident parcels plus\nthe element cytokines and debris (RK45):\n"
    "• recruitment ∝ D, split among parcels\n• polarization M0 ⇄ M1 ⇄ M2 (TNF-α, IL-10)\n"
    "• progenitor growth and differentiation\n• secretion and decay of c₁–c₄\n• phagocytosis of debris", num="1")
box(122, 60, 58, 36, "#e8f6f0", C2, "Transport",
    "Finite-volume diffusion of c₁–c₄\nbetween edge neighbours\n"
    "• mass-conserving, zero-flux boundaries\n• no decay in this step\n• debris does not move", num="2")
box(122, 17, 58, 36, "#fdeee6", C3, "Budding",
    "A secondary population > 1 cell-eq.\nbuds off as a new parcel in the same\nor the nearest free element\n"
    "• macrophage buds carry progenitors\n  (mixed-composition parcels)", num="3")
box(56, 17, 58, 36, "#fbf3dc", C4, "Migration",
    "Each parcel moves ≤ 1 element per hour\nby composition-weighted chemotaxis:\n"
    "• PMN, M0 → debris;  M1 → TNF-α\n• M2 → IL-10;  MSC → TGF-β, VEGF-A\n• softmax choice among free neighbours", num="4")
arrow((114.9, 78), (121.1, 78))
arrow((151, 59.2), (151, 53.8))
arrow((121.1, 35), (114.9, 35))
arrow((85, 53.8), (85, 59.2))
ax.add_patch(FancyBboxPatch((108.5, 53.6), 19, 6.0, boxstyle="round,pad=0.3,rounding_size=2", fc="white", ec=INK, lw=1.2))
ax.text(118, 56.6, "t → t + 1 h", ha="center", va="center", fontsize=8.5, color=INK, fontweight="bold")
arrow((47.0, 76), (55.2, 76), color=MUTED, lw=1.4)
ax.text(51.0, 78.5, "state", ha="center", fontsize=6.8, color=MUTED)
arrow((55.2, 42), (47.0, 42), color=MUTED, lw=1.4)
ax.text(51.0, 39.5, "new\npositions\n& parcels", ha="center", va="top", fontsize=6.0, color=MUTED, linespacing=1.15)

# ---------------------------------------------------------------- bottom: modules outside the loop
box(56, 1, 58, 12, "#f3f3f1", "#9a9a95", "Endothelial and oxygen modules",
    "implemented; disabled in this study", ls="--")
box(122, 1, 58, 12, "#f3f3f1", "#9a9a95", "Output to later healing",
    "fields and cells → mechanoregulatory model", ls="--")

fig.savefig(OUT / "fig6_schematic.png", dpi=600, bbox_inches="tight", facecolor="white")
fig.savefig(OUT / "fig6_schematic.pdf", bbox_inches="tight", facecolor="white")
print("written fig6_schematic")
