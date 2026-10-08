"""
Regression harness: run one 120-h simulation and print a JSON record of every
quantity the manuscript reports (calibration ratios, cytokine fold-changes,
zone composition, MSC by parcel type, parcel counts).

Usage:
    python tools/harness.py --params data/input_params.json --seed 1000 [--scripts scripts]
                            [--set K_lm=10.06 --set debris_coeff=0.5]

Seeding: if DomainModel accepts `seed=` (v2) it is passed directly. Otherwise
(v1.0.x code) every random.Random instance created during the run is forced to
`seed`, the same scheme used for the original calibration.
"""
import argparse
import json
import random
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
TARGET_HOURS = (72, 120)
POPULATIONS = {"PMN": 0, "M0": 1, "M1": 2, "M2": 3, "MSC": 7}
CYTOKINES = {"c1": 4, "c2": 5, "c3": 6, "c4": 8}


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--params", default=str(REPO / "data" / "input_params.json"))
    ap.add_argument("--mesh", default=str(REPO / "data" / "node_elements.txt"))
    ap.add_argument("--scripts", default=str(REPO / "scripts"))
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--set", action="append", default=[], help="override a parameter: name=value")
    ap.add_argument("--label", default="")
    ap.add_argument("--max_parcels", type=int, default=0, help="abort (screening) if parcels exceed this; 0 = off")
    return ap.parse_args()


def build_model(DomainModel, nodes, elements, params, seed):
    try:
        return DomainModel(nodes, elements, params, seed=seed)
    except TypeError:
        np.random.seed(seed)
        random.seed(seed)
        original = random.Random

        class _Seeded(original):
            def __init__(self, x=None):
                super().__init__(seed)

        random.Random = _Seeded
        try:
            return DomainModel(nodes, elements, params)
        finally:
            random.Random = original


def snapshot(model, ElementAgent, zone_elements):
    parcels = [a for a in model.schedule.agents if isinstance(a, ElementAgent)]
    by_element = {}
    for a in parcels:
        by_element.setdefault(a.element_id, []).append(a)

    totals = {k: float(sum(a.state[i] for a in parcels)) for k, i in {**POPULATIONS, **CYTOKINES}.items()}
    if hasattr(model, "cytokine_amounts"):   # v2: cytokines are element fields, not per-parcel copies
        totals.update(model.cytokine_amounts())
    totals["n_clipped"] = getattr(model, "n_clipped", 0)
    if hasattr(model, "cytokine_production"):
        totals.update({f"prod_{k}": v for k, v in model.cytokine_production().items()})
    totals["debris"] = float(sum(model.debris_field.values()))

    zones = {}
    for zone, eids in zone_elements.items():
        z = {k: 0.0 for k in POPULATIONS}
        for eid in eids:
            for a in by_element.get(eid, []):
                for k, i in POPULATIONS.items():
                    z[k] += a.state[i]
        z["M2_frac_M0M2"] = z["M2"] / (z["M0"] + z["M2"]) if z["M0"] + z["M2"] > 0 else None
        z["M2_frac_macro"] = z["M2"] / (z["M0"] + z["M1"] + z["M2"]) if z["M0"] + z["M1"] + z["M2"] > 0 else None
        zones[zone] = z

    by_type = {}
    for a in parcels:
        t = by_type.setdefault(a.agent_type, {"n": 0, **{k: 0.0 for k in POPULATIONS}})
        t["n"] += 1
        for k, i in POPULATIONS.items():
            t[k] += float(a.state[i])

    occupancy = [len(v) for v in model.element_agents.values()]
    return {"totals": totals, "zones": zones, "by_parcel_type": by_type,
            "n_parcels": len(parcels), "max_occupancy_list": max(occupancy)}


def ratio(a, b):
    return a / b if b else None


def main():
    args = parse_args()
    sys.path.insert(0, str(Path(args.scripts).resolve()))
    from utils import parse_nodes_elements
    from domain_model import DomainModel
    from element_agent_optimized import ElementAgent
    from zones import element_zones

    params = json.load(open(args.params))
    for item in args.set:
        name, value = item.split("=", 1)
        params[name] = float(value)

    nodes, elements = parse_nodes_elements(args.mesh)
    t0 = time.time()
    model = build_model(DomainModel, nodes, elements, params, args.seed)
    model.enable_EC = False
    zone_elements = element_zones(model.element_centroids)

    snaps = {"0": snapshot(model, ElementAgent, zone_elements)}
    while model.time < max(TARGET_HOURS):
        model.step()
        if args.max_parcels and len(model.schedule.agents) > args.max_parcels:
            print(json.dumps({"label": args.label, "seed": args.seed, "aborted": True,
                              "hour": model.time, "n_parcels": len(model.schedule.agents)}))
            return
        if model.time in TARGET_HOURS:
            snaps[str(model.time)] = snapshot(model, ElementAgent, zone_elements)

    s72, s120 = snaps["72"]["totals"], snaps["120"]["totals"]
    record = {
        "label": args.label, "seed": args.seed, "overrides": args.set,
        "elapsed_s": round(time.time() - t0, 1),
        "calibration": {
            "M1_M0_72": ratio(s72["M1"], s72["M0"]), "M2_M0_72": ratio(s72["M2"], s72["M0"]),
            "M1_M0_120": ratio(s120["M1"], s120["M0"]), "M2_M0_120": ratio(s120["M2"], s120["M0"]),
        },
        "cytokine_fold": {c: ratio(s120[c], s72[c]) for c in CYTOKINES},
        "production_fold": ({c: ratio(s120[f"prod_{c}"], s72[f"prod_{c}"]) for c in CYTOKINES}
                            if "prod_c1" in s72 else None),
        "snapshots": snaps,
    }
    print(json.dumps(record))


if __name__ == "__main__":
    main()
