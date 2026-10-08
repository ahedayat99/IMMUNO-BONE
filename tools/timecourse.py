"""
Full time-course record for one seed (figures and reviewer analyses):
  hourly : domain totals (PMN, M0, M1, M2, MSC), cytokine amounts and production rates, debris, parcels
  zones  : every 12 h, M0/M1/M2/PMN/MSC totals in the fracture-centre, periosteal and whole-callus zones
  every 24 h : parcel composition (M0, M1, M2 per macrophage-holding parcel) and
               capacity ratios per occupied element (PMN, macrophage, MSC)
  snapshots  : per-element cells, cytokines and debris at --snap hours (for spatial figures)
Usage: python tools/timecourse.py --params <json> --seed 11000 [--snap 0,24,72,120] > out.json
"""
import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from utils import parse_nodes_elements            # noqa: E402
from domain_model import DomainModel              # noqa: E402
from element_agent_optimized import ElementAgent  # noqa: E402
from zones import element_zones                   # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--params", required=True)
ap.add_argument("--seed", type=int, required=True)
ap.add_argument("--snap", default="")
ap.add_argument("--set", action="append", default=[])
args = ap.parse_args()
p = json.load(open(args.params))
for item in args.set:
    k, v = item.split("=", 1)
    p[k] = float(v)
snap_hours = {int(h) for h in args.snap.split(",") if h}
nodes, elements = parse_nodes_elements(REPO / "data" / "node_elements.txt")
model = DomainModel(nodes, elements, p, seed=args.seed)
model.enable_EC = False
zone_of = {}
for z, eids in element_zones(model.element_centroids).items():
    for e in eids:
        zone_of.setdefault(e, []).append(z)
POP = {"PMN": 0, "M0": 1, "M1": 2, "M2": 3, "MSC": 7}


def parcels():
    return [a for a in model.schedule.agents if isinstance(a, ElementAgent)]


def snapshot():
    cells = {}
    for a in parcels():
        c = cells.setdefault(a.element_id, [0.0] * 5)
        for i, k in enumerate(POP.values()):
            c[i] += float(a.state[k])
    return {str(e): {"x": model.element_centroids[e][0], "y": model.element_centroids[e][1],
                     "cells": cells.get(e, [0.0] * 5), "debris": model.debris_field[e],
                     "c": [model.cytokine_fields[k][e] for k in ("c1", "c2", "c3", "c4")]}
            for e in model.elements}


out = {"seed": args.seed, "hourly": [], "zones": {}, "composition": {}, "capacity": {}, "snapshots": {}}


def record(t):
    P = parcels()
    row = {"t": t, "parcels": len(P), "debris": float(sum(model.debris_field.values()))}
    for k, i in POP.items():
        row[k] = float(sum(a.state[i] for a in P))
    row.update({f"amt_{k}": v for k, v in model.cytokine_amounts().items()})
    row.update({f"prod_{k}": v for k, v in model.cytokine_production().items()})
    out["hourly"].append(row)
    if t % 12 == 0:
        acc = {}
        for a in P:
            for z in zone_of.get(a.element_id, []):
                v = acc.setdefault(z, [0.0] * 5)
                for j, i in enumerate(POP.values()):
                    v[j] += float(a.state[i])
        out["zones"][t] = acc
    if t % 24 == 0:
        out["composition"][t] = [[float(a.state[1]), float(a.state[2]), float(a.state[3])]
                                 for a in P if a.state[1] + a.state[2] + a.state[3] > 1e-3]
        tot = {}
        for a in P:
            v = tot.setdefault(a.element_id, [0.0, 0.0, 0.0])
            v[0] += a.state[0]
            v[1] += a.state[1] + a.state[2] + a.state[3]
            v[2] += a.state[7]
        out["capacity"][t] = [[v[0] / (p["PMN_max"] * model.area_ratio[e]),
                               v[1] / (p["M_max"] * model.area_ratio[e]),
                               v[2] / (p["K_lm"] * model.area_ratio[e])] for e, v in tot.items()]
    if t in snap_hours:
        out["snapshots"][t] = snapshot()


record(0)
while model.time < 120:
    model.step()
    record(model.time)
out["n_clipped"] = model.n_clipped
print(json.dumps(out))
