"""
Distribution of macrophage phenotype composition across parcels (analysis for Reviewer Comment 2).
For every parcel holding macrophages, records M0, M1, M2 amounts at the requested hours.
Usage: python tools/parcel_composition.py --params <json> --seed 11000 > out.json
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
ap.add_argument("--hours", default="24,48,72,96,120")
args = ap.parse_args()
hours = {int(h) for h in args.hours.split(",")}
nodes, elements = parse_nodes_elements(REPO / "data" / "node_elements.txt")
model = DomainModel(nodes, elements, json.load(open(args.params)), seed=args.seed)
model.enable_EC = False
out = {"seed": args.seed, "t": {}}
while model.time < max(hours):
    model.step()
    if model.time in hours:
        out["t"][model.time] = [[a.agent_type, float(a.state[1]), float(a.state[2]), float(a.state[3]),
                                 float(a.state[0]), float(a.state[7])]
                                for a in model.schedule.agents if isinstance(a, ElementAgent)]
print(json.dumps(out))
