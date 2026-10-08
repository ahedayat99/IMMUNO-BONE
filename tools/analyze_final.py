"""
Summaries of the model-of-record suite (out/final/timecourse): perturbation effect sizes at 120 h,
baseline time courses, zone M2 enrichment, parcel composition and capacity exceedance.
Writes out/final/summary.json and prints the key tables.
"""
import glob
import json
import re
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parents[1] / "out" / "final"
runs = {}
for f in sorted(glob.glob(str(OUT / "timecourse" / "*.json"))):
    m = re.match(r"(.+)_(\d+)\.json", Path(f).name)
    runs.setdefault(m.group(1), []).append(json.load(open(f)))

KEYS = ["PMN", "M0", "M1", "M2", "MSC", "amt_c1", "amt_c2", "amt_c3", "amt_c4", "debris"]


def at(run, t):
    return next(r for r in run["hourly"] if r["t"] == t)


summary = {"n_seeds": {k: len(v) for k, v in runs.items()}}

# ---- perturbation effect sizes at 120 h (% change of seed-mean vs baseline, same seeds)
base = {k: np.mean([at(r, 120)[k] for r in runs["baseline"]]) for k in KEYS}
eff = {}
for name, R in runs.items():
    if name == "baseline":
        continue
    m = {k: np.mean([at(r, 120)[k] for r in R]) for k in KEYS}
    eff[name] = {k: 100 * (m[k] - base[k]) / base[k] for k in KEYS}
summary["effect_120h_percent"] = eff
print(f"{'scenario':14s}" + "".join(f"{k.replace('amt_', ''):>9s}" for k in KEYS) + "   n")
for name, e in eff.items():
    print(f"{name:14s}" + "".join(f"{e[k]:+9.1f}" for k in KEYS) + f"  {len(runs[name])}")

# ---- baseline time courses (mean, SD per hour)
B = runs["baseline"]
tc = {k: {"mean": [float(np.mean([r["hourly"][t][k] for r in B])) for t in range(121)],
          "sd": [float(np.std([r["hourly"][t][k] for r in B])) for t in range(121)]}
      for k in KEYS + ["prod_c1", "prod_c2", "prod_c3", "prod_c4", "parcels"]}
summary["baseline_timecourse"] = tc
peaks = {k: int(np.argmax(tc[k]["mean"])) for k in ("PMN", "M0", "M1", "M2", "MSC")}
summary["peak_hour"] = peaks
print("\npeak hour of each population (baseline):", peaks)
d = tc["debris"]["mean"]
summary["debris_remaining_percent"] = {t: 100 * d[t] / d[0] for t in (24, 48, 72, 96, 120)}
print("debris remaining (% of initial):", {t: round(v, 1) for t, v in summary["debris_remaining_percent"].items()})

# ---- zone M2 enrichment (pooled sums across seeds), M2/(M0+M2)
zr = {}
for t in sorted({int(k) for k in B[0]["zones"]}):
    def frac(z):
        s = np.sum([r["zones"][str(t)].get(z, [0] * 5) for r in B], axis=0)
        return s[3] / (s[1] + s[3]) if s[1] + s[3] > 0 else np.nan
    c, p = frac("fracture_centre"), frac("periosteal")
    per_seed = []
    for r in B:
        zc, zp = r["zones"][str(t)].get("fracture_centre"), r["zones"][str(t)].get("periosteal")
        if zc and zp and zc[1] + zc[3] > 0 and zp[1] + zp[3] > 0:
            per_seed.append((zp[3] / (zp[1] + zp[3])) / (zc[3] / (zc[1] + zc[3])))
    zr[t] = {"centre": c, "periosteal": p, "ratio_pooled": p / c if c else np.nan,
             "seeds_ratio_gt1": int(sum(x > 1 for x in per_seed)), "n": len(per_seed)}
summary["zones"] = zr
print("\nperiosteal/centre M2 fraction (pooled):", {t: round(v["ratio_pooled"], 2) for t, v in zr.items()})
print("seeds with ratio > 1:", {t: f"{v['seeds_ratio_gt1']}/{v['n']}" for t, v in zr.items()})

# ---- parcel composition and capacity
comp, cap = {}, {}
for t in ("24", "48", "72", "96", "120"):
    M = np.array([row for r in B for row in r["composition"][t]])
    fr = M / M.sum(1, keepdims=True)
    comp[t] = {"parcels_per_seed": len(M) / len(B), "mixed_fraction": float(np.mean(fr.max(1) <= 0.95)),
               "M2_share_p10_p50_p90": np.percentile(fr[:, 2], [10, 50, 90]).tolist(),
               "M2_share_hist10": np.histogram(fr[:, 2], bins=10, range=(0, 1))[0].tolist()}
    C = np.array([row for r in B for row in r["capacity"][t]])
    cap[t] = {"occupied_per_seed": len(C) / len(B),
              "PMN_over": float(np.mean(C[:, 0] > 1)), "PMN_max_ratio": float(C[:, 0].max()),
              "mac_over": float(np.mean(C[:, 1] > 1)), "mac_max_ratio": float(C[:, 1].max()),
              "MSC_over": float(np.mean(C[:, 2] > 1)), "MSC_max_ratio": float(C[:, 2].max())}
summary["composition"], summary["capacity"] = comp, cap
print("\nmixed parcels:", {t: round(v["mixed_fraction"], 3) for t, v in comp.items()},
      " M2 share p10/50/90 at 120 h:", np.round(comp["120"]["M2_share_p10_p50_p90"], 3))
print("capacity exceedance:", {t: (round(v["PMN_over"], 3), round(v["mac_over"], 4), round(v["MSC_over"], 4)) for t, v in cap.items()})
mac = [sum(at(r, 120)[k] for k in ("M0", "M1", "M2")) / sum(at(r, 72)[k] for k in ("M0", "M1", "M2")) for r in B]
summary["unfitted_macrophage_fold"] = [float(np.mean(mac)), float(np.std(mac))]
print("unfitted total macrophage fold D5/D3: %.2f ± %.2f" % tuple(summary["unfitted_macrophage_fold"]))
json.dump(summary, open(OUT / "summary.json", "w"), indent=1, default=float)
