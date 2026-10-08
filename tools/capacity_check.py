"""
Comment 6: how often, and by how much, element totals exceed their lineage capacities
(PMN_max, M_max for M0+M1+M2, K_lm for MSC; each scaled by the element's area ratio a_e).
Usage: python tools/capacity_check.py --params <json> --seed 11000 > out.json
"""
import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from utils import parse_nodes_elements          # noqa: E402
from domain_model import DomainModel            # noqa: E402
from element_agent_optimized import ElementAgent  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--params", required=True)
ap.add_argument("--seed", type=int, required=True)
args = ap.parse_args()
p = json.load(open(args.params))
nodes, elements = parse_nodes_elements(REPO / "data" / "node_elements.txt")
model = DomainModel(nodes, elements, p, seed=args.seed)
model.enable_EC = False
out = {"seed": args.seed, "t": {}}
while model.time < 120:
    model.step()
    if model.time % 24 == 0:
        tot = {}
        for a in model.schedule.agents:
            if isinstance(a, ElementAgent):
                t = tot.setdefault(a.element_id, [0.0, 0.0, 0.0])
                t[0] += a.state[0]
                t[1] += a.state[1] + a.state[2] + a.state[3]
                t[2] += a.state[7]
        rows = []
        for eid, (N, M, C) in tot.items():
            ae = model.area_ratio[eid]
            rows.append([N / (p["PMN_max"] * ae), M / (p["M_max"] * ae), C / (p["K_lm"] * ae), N, M, C])
        out["t"][model.time] = rows
print(json.dumps(out))
