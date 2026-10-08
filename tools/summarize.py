"""
Summarize one or more harness JSONL files side by side (mean ± SD over seeds).

Usage:
    python tools/summarize.py tests/golden/v1_baseline.jsonl out/v2_step.jsonl
"""
import json
import sys

import numpy as np

# Experimental targets (in-vivo_outputs/IF.xlsx "used these (final)"; dPCR fold-changes as used in calibration)
TARGETS = {
    "M1_M0_72": 1.620, "M2_M0_72": 0.145, "M1_M0_120": 1.667, "M2_M0_120": 1.367,
    "c1": 12.85, "c2": 6.20, "c3": 8.97, "c4": 5.54,
}


def load(path):
    return [json.loads(line) for line in open(path) if line.strip()]


def rows(records):
    out = {}
    for k in ("M1_M0_72", "M2_M0_72", "M1_M0_120", "M2_M0_120"):
        out[k] = [r["calibration"][k] for r in records]
    for c in ("c1", "c2", "c3", "c4"):
        out[f"{c} fold"] = [r["cytokine_fold"][c] for r in records]
    for t in ("72", "120"):
        for pop in ("PMN", "M0", "M1", "M2", "MSC"):
            out[f"{pop} total {t}h"] = [r["snapshots"][t]["totals"][pop] for r in records]
        z = [r["snapshots"][t]["zones"] for r in records]
        out[f"centre M2/(M0+M2) {t}h"] = [s["fracture_centre"]["M2_frac_M0M2"] for s in z]
        out[f"perio M2/(M0+M2) {t}h"] = [s["periosteal"]["M2_frac_M0M2"] for s in z]
        out[f"perio/centre {t}h"] = [
            (s["periosteal"]["M2_frac_M0M2"] / s["fracture_centre"]["M2_frac_M0M2"])
            if s["periosteal"]["M2_frac_M0M2"] and s["fracture_centre"]["M2_frac_M0M2"] else None
            for s in z]
        out[f"parcels {t}h"] = [r["snapshots"][t]["n_parcels"] for r in records]
    return out


def fmt(vals):
    v = np.array([x for x in vals if x is not None], dtype=float)
    if v.size == 0:
        return "n/a"
    return f"{v.mean():.3f} ± {v.std():.3f}"


def main():
    files = sys.argv[1:]
    tables = [rows(load(f)) for f in files]
    head = f"{'quantity':26s}{'target':>8s}" + "".join(f"{f.split('/')[-1][:24]:>26s}" for f in files)
    print(head)
    for key in tables[0]:
        tk = key.replace(" fold", "")
        tgt = f"{TARGETS[tk]:8.3f}" if tk in TARGETS else f"{'':8s}"
        print(f"{key:26s}{tgt}" + "".join(f"{fmt(t[key]):>26s}" for t in tables))


if __name__ == "__main__":
    main()
