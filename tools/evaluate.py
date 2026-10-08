"""
Shared evaluation for sensitivity analysis and calibration.

evaluate(params_file, overrides, seeds) runs tools/harness.py once per seed (separate processes)
and returns seed-averaged outputs:
  ratios     : M1/M0 and M2/M0 at 72 h and 120 h (immunofluorescence targets)
  prod_fold  : Day 5 / Day 3 cytokine production-rate fold (compared with dPCR transcript folds)
  conc_fold  : Day 5 / Day 3 cytokine amount fold (reported for comparison with v1.0.1)
  perio_ctr  : periosteal / fracture-centre M2 fraction at 72 h and 120 h
Runs whose parcel count exceeds max_parcels are aborted and scored with PENALTY.
"""
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
TARGETS_M = {"M1_M0_72": 1.620, "M2_M0_72": 0.145, "M1_M0_120": 1.667, "M2_M0_120": 1.367}
TARGETS_C = {"c1": 12.85, "c2": 6.20, "c3": 8.97, "c4": 5.54}
PENALTY = 1e3
OUTPUTS = list(TARGETS_M) + [f"prod_{c}" for c in TARGETS_C]   # the 8 experimentally anchored outputs


def run_one(params_file, overrides, seed, max_parcels=3000, timeout=3600):
    cmd = [sys.executable, str(REPO / "tools" / "harness.py"), "--seed", str(seed),
           "--params", str(params_file), "--max_parcels", str(max_parcels)]
    for k, v in overrides.items():
        cmd += ["--set", f"{k}={float(v):.10g}"]
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, env=env)
    except subprocess.TimeoutExpired:
        return {"failed": "timeout"}
    if r.returncode != 0:
        return {"failed": r.stderr[-300:]}
    rec = json.loads(r.stdout.strip().splitlines()[-1])
    if rec.get("aborted"):
        return {"failed": f"aborted: {rec['n_parcels']} parcels at {rec['hour']} h"}
    z = {t: rec["snapshots"][t]["zones"] for t in ("72", "120")}

    def pc(t):
        c, p = z[t]["fracture_centre"]["M2_frac_M0M2"], z[t]["periosteal"]["M2_frac_M0M2"]
        return p / c if c and p is not None else None

    out = dict(rec["calibration"])
    out.update({f"prod_{c}": rec["production_fold"][c] for c in TARGETS_C})
    out.update({f"conc_{c}": rec["cytokine_fold"][c] for c in TARGETS_C})
    out.update({"perio_ctr_72": pc("72"), "perio_ctr_120": pc("120"),
                "parcels_120": rec["snapshots"]["120"]["n_parcels"], "elapsed_s": rec["elapsed_s"]})
    return out


def average(records):
    ok = [r for r in records if "failed" not in r]
    if len(ok) < len(records):
        return {"failed": [r["failed"] for r in records if "failed" in r][0]}
    keys = ok[0].keys()
    return {k: float(np.mean([r[k] for r in ok if r[k] is not None])) if any(r[k] is not None for r in ok) else None
            for k in keys}


def loss_macrophage(out):
    """Calibration objective: sum of squared relative errors of the four macrophage ratios."""
    if "failed" in out:
        return PENALTY
    return float(sum(((out[k] - t) / t) ** 2 for k, t in TARGETS_M.items()))


def cytokine_directions(out):
    """Validation criterion: number of cytokines whose production rises from Day 3 to Day 5 (dPCR: all four rise)."""
    return sum(out[f"prod_{c}"] > 1 for c in TARGETS_C) if "failed" not in out else 0
