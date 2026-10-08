"""
Build the sensitivity-analysis base point: option parameters + the best pilot sample (macrophage loss).
Usage: python tools/make_sa_base.py --option A --pilot out/option_A/pilot_lhs.jsonl
"""
import argparse
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ap = argparse.ArgumentParser()
ap.add_argument("--option", choices=["A", "B"], required=True)
ap.add_argument("--pilot", required=True)
args = ap.parse_args()
R = [json.loads(l) for l in open(args.pilot) if l.strip()]
best = min((r for r in R if "calibration" in r), key=lambda r: r["loss_M"])
base = json.load(open(REPO / "params" / f"option_{args.option}" / "params.json"))
base.update(best["params"])
json.dump(base, open(REPO / "params" / f"option_{args.option}" / "sa_base.json", "w"), indent=4)
print(f"sa_base.json from pilot sample {best['i']} (loss_M {best['loss_M']:.3f}); calibration {best['calibration']}")
