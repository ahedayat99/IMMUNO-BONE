"""
Comment 9: convergence of the finite-volume transport operator under mesh refinement.
Each eight-node quadrilateral is split into four (corner, mid-side and centre nodes). A Gaussian is diffused
on the native and refined meshes and compared with the exact free-space solution in an interior window.
"""
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from utils import parse_nodes_elements  # noqa: E402
from mesh import MeshGeometry           # noqa: E402


def refine(nodes, elements):
    nodes = dict(nodes)
    nid = max(nodes) + 1
    new = {}
    eid = 1
    for _, n in elements.items():
        c1, c2, c3, c4, m12, m23, m34, m41 = n
        xc = tuple(np.mean([nodes[k] for k in n[:4]], axis=0))
        nodes[nid] = xc
        cen = nid
        nid += 1
        for quad in ([c1, m12, cen, m41], [m12, c2, m23, cen], [cen, m23, c3, m34], [m41, cen, m34, c4]):
            new[eid] = quad                       # four-node quads (corners only)
            eid += 1
    return nodes, new


def gaussian(xy, x0, D, t):
    r2 = np.sum((xy - x0) ** 2, axis=1)
    return np.exp(-r2 / (4 * D * t)) / (4 * np.pi * D * t)


nodes, elements = parse_nodes_elements(REPO / "data" / "node_elements.txt")
meshes = {"native": MeshGeometry(nodes, elements), "refined x4": MeshGeometry(*refine(nodes, elements))}
D, t0, T = 0.108, 0.5, 1.0
for x0 in (np.array([4.5, 5.0]), np.array([3.5, -5.0])):
    print(f"source at {x0.tolist()}:")
    for name, g in meshes.items():
        c = gaussian(g.centroid, x0, D, t0)
        c = g.diffuse(c, D, T)
        exact = gaussian(g.centroid, x0, D, t0 + T)
        w = (np.abs(g.centroid[:, 0] - x0[0]) < 1.5) & (np.abs(g.centroid[:, 1] - x0[1]) < 1.5)
        err = np.sqrt(np.sum(g.area[w] * (c[w] - exact[w]) ** 2) / np.sum(g.area[w] * exact[w] ** 2))
        h = np.sqrt(g.mean_area)
        print(f"  {name:11s} elements {len(g.ids):6d}  mean size {h:.3f} mm  relative L2 error {100*err:.2f}%")
