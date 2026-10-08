"""
Summarize an LHS pilot: best samples by macrophage-only loss and by the 8-target loss,
what those samples predict for cytokines and zones, and where good samples sit in parameter space.

Usage: python tools/analyze_pilot.py out/option_B/pilot_lhs.jsonl
"""
import json
import sys

import numpy as np

TARGETS_M = {"M1_M0_72": 1.620, "M2_M0_72": 0.145, "M1_M0_120": 1.667, "M2_M0_120": 1.367}
TARGETS_C = {"c1": 12.85, "c2": 6.20, "c3": 8.97, "c4": 5.54}

R = [json.loads(l) for l in open(sys.argv[1]) if l.strip()]
ok = [r for r in R if "error" not in r and not r.get("aborted") and np.isfinite(r["loss_M"])]
print(f"samples {len(R)}, finished {len(ok)}, errors {len(R) - len(ok)}")


def mape(r, keys, targets, sec):
    return 100 * np.mean([abs(r[sec][k] - targets[k]) / targets[k] for k in keys])


def perio_ratio(r, t):
    z = r[f"zones_{t}"]
    c, p = z["fracture_centre"]["M2_frac_M0M2"], z["periosteal"]["M2_frac_M0M2"]
    return p / c if c and p is not None else float("nan")


for loss in ("loss_M", "loss_all"):
    best = sorted(ok, key=lambda r: r[loss])[:8]
    print(f"\n=== best 8 by {loss} ===")
    print(f"{'loss':>8} {'M1/M0 72':>9} {'M2/M0 72':>9} {'M1/M0 120':>10} {'M2/M0 120':>10} | "
          f"{'TNF':>6} {'IL10':>6} {'TGF':>6} {'VEGF':>6} | {'MAPE_M%':>7} {'perio/ctr 72':>12} {'120':>5} {'parcels':>7}")
    for r in best:
        c, f = r["calibration"], r["cytokine_fold"]
        print(f"{r[loss]:8.3f} {c['M1_M0_72']:9.3f} {c['M2_M0_72']:9.3f} {c['M1_M0_120']:10.3f} {c['M2_M0_120']:10.3f} | "
              f"{f['c1']:6.2f} {f['c2']:6.2f} {f['c3']:6.2f} {f['c4']:6.2f} | "
              f"{mape(r, TARGETS_M, TARGETS_M, 'calibration'):7.1f} {perio_ratio(r, '72'):12.2f} {perio_ratio(r, '120'):5.2f} {r['parcels_120']:7d}")
    print("targets  1.620     0.145      1.667      1.367     | 12.85   6.20   8.97   5.54")

print("\n=== fraction of samples with correct cytokine directions (all four folds > 1) ===")
dirs = [all(r["cytokine_fold"][k] > 1 for k in TARGETS_C) for r in ok]
print(f"{np.mean(dirs):.2f}")

good = sorted(ok, key=lambda r: r["loss_M"])[: max(10, len(ok) // 20)]
print(f"\n=== parameter ranges of the best {len(good)} samples by loss_M (log10) ===")
for k in ok[0]["params"]:
    v = np.log10([r["params"][k] for r in good])
    allv = np.log10([r["params"][k] for r in ok])
    print(f"{k:12s} best: {v.min():6.2f} .. {v.max():6.2f} (median {np.median(v):6.2f})   sampled: {allv.min():6.2f} .. {allv.max():6.2f}")
