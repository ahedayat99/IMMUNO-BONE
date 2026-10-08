"""Maximum relative difference between two harness JSONL files, seed by seed."""
import json
import sys

a = {r["seed"]: r for r in map(json.loads, open(sys.argv[1]))}
b = {r["seed"]: r for r in map(json.loads, open(sys.argv[2]))}
worst, where = 0.0, None
for seed in a:
    x, y = a[seed], b[seed]
    pairs = [(f"{sec}.{k}", x[sec][k], y[sec][k]) for sec in ("calibration", "cytokine_fold") for k in x[sec]]
    for t in ("72", "120"):
        pairs += [(f"{t}h.{k}", v, y["snapshots"][t]["totals"][k]) for k, v in x["snapshots"][t]["totals"].items()]
    for name, u, v in pairs:
        d = abs(u - v) / (abs(u) + 1e-12)
        if d > worst:
            worst, where = d, (seed, name)
print(f"max relative difference: {worst:.2e} at {where}")
