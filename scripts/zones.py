"""
Anatomical zones of the callus mesh, matched to the histological regions of interest.

Coordinates are in mm, origin at the centre of the osteotomy gap.
x is perpendicular to the fixation plate (layers along +x: plate | cis cortex |
medullary canal | trans cortex | callus); y runs along the bone axis.
Only the osteotomy gap and the periosteal callus in front of the trans cortex are meshed.
"""

GAP_OUTER_X = 1.88      # boundary between the osteotomy gap and the callus (mm)
GAP_HALF_HEIGHT = 1.0   # the 2-mm osteotomy spans |y| <= 1 mm
PERIOSTEAL_Y = (1.0, 8.0)  # histological periosteal ROI (QuPath): from the gap edge along the callus to |y| = 8 mm
PERIOSTEAL_DEPTH = 2.0     # ... and 2 mm deep into the callus from the far cortex (x = 1.88-3.88 mm)
# Submitted-revision definition, kept for reference: x >= 1.88, 2 <= |y| <= 3 (full callus depth).

ZONES = {
    "fracture_centre": lambda x, y: x < GAP_OUTER_X and abs(y) <= GAP_HALF_HEIGHT,
    "periosteal": lambda x, y: (GAP_OUTER_X <= x <= GAP_OUTER_X + PERIOSTEAL_DEPTH
                               and PERIOSTEAL_Y[0] <= abs(y) <= PERIOSTEAL_Y[1]),
    "whole_callus": lambda x, y: True,
}


def element_zones(element_centroids):
    """Map each zone name to the list of element ids whose centroid lies in it."""
    return {
        name: [eid for eid, (x, y) in element_centroids.items() if inside(x, y)]
        for name, inside in ZONES.items()
    }
