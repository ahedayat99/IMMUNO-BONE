"""
Numerical verification of the molecular-field operators (Supplementary Note S4).
  1. Mass conservation of finite-volume transport on the callus mesh.
  2. Cytokine decay in the reaction step (no cells, no transport) against exp(-d t).
  3. Convergence of transport on regular quadrilateral meshes (diffusing Gaussian, exact solution).
  4. Non-orthogonality of the callus mesh (angle between the centroid link and the shared-edge normal).
Usage: python tools/transport_verification.py   (prints a JSON summary)
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from utils import parse_nodes_elements    # noqa: E402
from mesh import MeshGeometry             # noqa: E402
from reaction import react_element        # noqa: E402

p = json.load(open(REPO / "params" / "model_of_record.json"))
nodes, elements = parse_nodes_elements(REPO / "data" / "node_elements.txt")
g = MeshGeometry(nodes, elements)
out = {}

# 1. mass conservation
rng = np.random.default_rng(0)
c0 = rng.random(len(g.ids))
c1 = g.diffuse(c0.copy(), p["diff_c1"], 1.0)
out["mass_relative_change"] = float(abs(np.sum(g.area * c1) - np.sum(g.area * c0)) / np.sum(g.area * c0))

# 2. decay in the reaction step (empty element, no debris)
f0 = np.array([0.0, 1.0, 1.0, 1.0, 1.0])
f1, _, _ = react_element(f0, np.zeros((0, 5)), p, 1.0)
exact = np.exp(-np.array([p["d_c1"], p["d_c2"], p["d_c3"], p["d_c4"]]))
out["decay_max_relative_error"] = float(np.max(np.abs(f1[1:] - exact) / exact))


# 3. regular meshes
def regular(h, L=3.0):
    n = int(round(2 * L / h))
    nd, el, k = {}, {}, 1
    idx = lambda i, j: i * (n + 1) + j + 1
    for i in range(n + 1):
        for j in range(n + 1):
            nd[idx(i, j)] = (-L + i * h, -L + j * h)
    for i in range(n):
        for j in range(n):
            el[k] = [idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)]
            k += 1
    return MeshGeometry(nd, el)


def gaussian(xy, D, t):
    return np.exp(-np.sum(xy ** 2, axis=1) / (4 * D * t)) / (4 * np.pi * D * t)


D, t0, T = p["diff_c1"], 0.5, 1.0
conv = []
for h in (0.2, 0.1, 0.05):
    m = regular(h)
    c = m.diffuse(gaussian(m.centroid, D, t0), D, T)
    ex = gaussian(m.centroid, D, t0 + T)
    w = np.all(np.abs(m.centroid) < 1.5, axis=1)
    err = math.sqrt(np.sum(m.area[w] * (c[w] - ex[w]) ** 2) / np.sum(m.area[w] * ex[w] ** 2))
    conv.append({"h_mm": h, "elements": len(m.ids), "rel_L2_error_percent": 100 * err})
out["regular_mesh"] = conv
out["observed_order"] = [math.log(conv[i]["rel_L2_error_percent"] / conv[i + 1]["rel_L2_error_percent"], 2)
                         for i in range(len(conv) - 1)]

# 4. non-orthogonality of the callus mesh
edge = {}
for i, eid in enumerate(g.ids):
    c = elements[eid][:4]
    for k in range(4):
        a, b = c[k], c[(k + 1) % 4]
        edge.setdefault((min(a, b), max(a, b)), []).append(i)
ang = []
for (a, b), owners in edge.items():
    if len(owners) != 2:
        continue
    t = np.subtract(nodes[b], nodes[a])
    nrm = np.array([t[1], -t[0]]) / np.linalg.norm(t)
    d = g.centroid[owners[1]] - g.centroid[owners[0]]
    ang.append(math.degrees(math.acos(min(1.0, abs(np.dot(nrm, d)) / np.linalg.norm(d)))))
ang = np.array(ang)
out["non_orthogonality_deg"] = {"median": float(np.median(ang)), "p95": float(np.percentile(ang, 95)),
                                "max": float(ang.max()), "interior_edges": int(len(ang))}
print(json.dumps(out, indent=2))
