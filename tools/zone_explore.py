"""
Zone-definition sensitivity of the periosteal vs fracture-centre M2 enrichment (analysis only).
Records, every 12 h, M0/M1/M2 totals in several candidate zone definitions.
Usage: python tools/zone_explore.py --params params/selected_v2/params.json --seed 11000 > out.json
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

G = 1.88  # gap / callus boundary (mm)
ZONES = {
    "centre_gap":            lambda x, y: x < G and abs(y) <= 1.0,
    "centre_gap+front":      lambda x, y: abs(y) <= 1.0,
    "front_of_gap":          lambda x, y: x >= G and abs(y) <= 1.0,
    "perio_y2-3_full":       lambda x, y: x >= G and 2.0 <= abs(y) <= 3.0,    # current definition (author: 1-2 mm beyond gap)
    "perio_y2-3_band1mm":    lambda x, y: G <= x <= G + 1.0 and 2.0 <= abs(y) <= 3.0,
    "perio_y1-3_full":       lambda x, y: x >= G and 1.0 < abs(y) <= 3.0,
    "perio_y1-2_full":       lambda x, y: x >= G and 1.0 < abs(y) <= 2.0,
    "perio_y1-3_band1mm":    lambda x, y: G <= x <= G + 1.0 and 1.0 < abs(y) <= 3.0,
    "callus_outside_gap_y":  lambda x, y: x >= G and abs(y) > 1.0,
}

ap = argparse.ArgumentParser()
ap.add_argument("--params", required=True)
ap.add_argument("--seed", type=int, required=True)
args = ap.parse_args()
nodes, elements = parse_nodes_elements(REPO / "data" / "node_elements.txt")
model = DomainModel(nodes, elements, json.load(open(args.params)), seed=args.seed)
model.enable_EC = False
membership = {eid: [z for z, f in ZONES.items() if f(x, y)] for eid, (x, y) in model.element_centroids.items()}
out = {"seed": args.seed, "t": {}}
while model.time < 120:
    model.step()
    if model.time % 12 == 0:
        acc = {z: [0.0, 0.0, 0.0] for z in ZONES}
        for a in model.schedule.agents:
            if isinstance(a, ElementAgent):
                for z in membership[a.element_id]:
                    for k in range(3):
                        acc[z][k] += float(a.state[1 + k])
        out["t"][model.time] = acc
print(json.dumps(out))
