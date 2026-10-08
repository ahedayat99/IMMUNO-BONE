"""
Figure 6: keep the original graphic (figures_v5.pptx, slide 6, image28.png) and replace only its text.
The original is upscaled ×3, the old text blocks and the watermark are painted over with the white
background, and corrected text is drawn in Helvetica Neue at the same positions.
Usage: python tools/fig6_retext.py <original.png>  -> out/figures/fig6_architecture.png (+ .pdf)
"""
import sys
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

REPO = Path(__file__).resolve().parents[1]
src = Image.open(sys.argv[1]).convert("RGB")
K = 3                                                     # upscale factor
im = src.resize((src.width * K, src.height * K), Image.LANCZOS)
d = ImageDraw.Draw(im)
HN = "/System/Library/Fonts/HelveticaNeue.ttc"
BOLD, MED, REG = 1, 10, 0


def font(size, face):
    return ImageFont.truetype(HN, size * K, index=face)


def erase(x0, y0, x1, y1):
    d.rectangle([x0 * K, y0 * K, x1 * K, y1 * K], fill="white")


def text(cx, y, lines, size, face, color, gap=1.22):
    f = font(size, face)
    for i, ln in enumerate(lines):
        w = d.textlength(ln, font=f)
        d.text((cx * K - w / 2, (y + i * size * gap) * K), ln, font=f, fill=color)


TITLE, HEAD, BODY, CENTRE = "#4f5453", "#141414", "#2b2b2b", "#555b5a"

# ---- erase old text (coordinates in the original 1376 x 768 image)
erase(160, 20, 1220, 95)          # title
erase(555, 222, 825, 294)         # top node
erase(800, 438, 1045, 510)        # right node
erase(345, 438, 565, 510)         # left node
erase(555, 666, 825, 742)         # bottom node
erase(598, 364, 778, 488)         # centre circle text
erase(1255, 735, 1376, 766)       # watermark
# clear the whole centre circle and redraw its outline (measured: centre (687.5, 420.3), radius 118.9)
CX, CY, R = 687.5, 420.3, 118.9
d.ellipse([(CX - R - 3) * K, (CY - R - 3) * K, (CX + R + 3) * K, (CY + R + 3) * K], fill="white")
d.ellipse([(CX - R) * K, (CY - R) * K, (CX + R) * K, (CY + R) * K], outline=(98, 101, 105), width=int(1.6 * K))

# ---- new text
text(688, 30, ["Architecture: Parcel–Field Coupling in One 1-h Step"], 40, BOLD, TITLE)

text(688, 228, ["Population Parcel"], 19, BOLD, HEAD)
text(688, 254, ["Mixed PMN, M0, M1, M2 and MSC;", "element ODEs: recruitment, polarization"], 13.5, REG, BODY)

text(920, 443, ["Secretion"], 19, BOLD, HEAD)
text(920, 469, ["Parcels secrete cytokines into", "fields shared by the element"], 13.5, REG, BODY)

text(688, 673, ["Continuous Fields"], 19, BOLD, HEAD)
text(688, 699, ["Finite-volume diffusion of", "cytokines on the FE mesh"], 13.5, REG, BODY)

text(455, 443, ["Budding & Chemotaxis"], 19, BOLD, HEAD)
text(455, 469, ["New mixed parcels bud; fields", "steer composition-weighted migration"], 13.5, REG, BODY)

text(688, 368, ["Operator-split", "coupling (Δt = 1 h):", "parcels set local", "sources; fields set", "polarization and", "migration."], 15.5, BOLD, CENTRE, gap=1.18)

out = REPO / "out" / "figures" / "fig6_architecture"
im.save(f"{out}.png", dpi=(600, 600))
im.save(f"{out}.pdf", resolution=600)
print("written", out)
