"""
One-at-a-time (OAT) sensitivity analysis, as in Methods 5.10.1:
  S_ij(eps) = [g_j(theta_i(1+eps)) - g_j(theta_i(1-eps))] / (2 eps g_j(theta)),   eps in {0.05, 0.10, 0.15}
  S_ij = mean over eps.
Outputs g_j: the eight experimentally anchored quantities (four macrophage ratios, four cytokine
production folds). Each evaluation is the mean of the same seeds (common random numbers).
Selection rule: parameters appearing in the top 10 by |S_ij| for >= min_outputs of the 8 outputs.

Usage:
    python tools/oat_sa.py --base params/option_B/sa_base.json --out out/option_B/sa --workers 32
"""
import argparse
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

from evaluate import OUTPUTS, average, run_one

EPS = (0.05, 0.10, 0.15)
# Initial-condition, scenario and numerical keys are not model parameters of Table S1
EXCLUDE_PREFIX = ("init_", "chemotaxis_")
EXCLUDE = {"debris_coeff", "max_parcels_per_element"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True, help="full parameter JSON at which sensitivities are evaluated")
    ap.add_argument("--out", required=True)
    ap.add_argument("--seeds", default="1000,2000,3000")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--top", type=int, default=10)
    ap.add_argument("--min_outputs", type=int, default=3)
    ap.add_argument("--names", default="", help="comma-separated subset (testing only)")
    args = ap.parse_args()

    base = json.load(open(args.base))
    names = [k for k, v in base.items()
             if isinstance(v, (int, float)) and k not in EXCLUDE and not k.startswith(EXCLUDE_PREFIX) and v != 0]
    if args.names:
        names = [k for k in names if k in args.names.split(",")]
    seeds = [int(s) for s in args.seeds.split(",")]
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    raw_file = out / "raw_runs.jsonl"
    done = {}
    if raw_file.exists():
        for line in open(raw_file):
            r = json.loads(line)
            done[(r["param"], r["factor"], r["seed"])] = r["out"]

    jobs = [("__base__", 1.0, s) for s in seeds]
    jobs += [(k, 1 + sgn * e, s) for k in names for e in EPS for sgn in (1, -1) for s in seeds]
    todo = [j for j in jobs if j not in done]
    print(f"{len(names)} parameters, {len(jobs)} runs, {len(todo)} to do")

    def run(job):
        k, f, s = job
        ov = {} if k == "__base__" else {k: base[k] * f}
        return job, run_one(args.base, ov, s)

    with ThreadPoolExecutor(max_workers=args.workers) as pool, open(raw_file, "a") as fh:
        for (k, f, s), res in pool.map(run, todo):
            done[(k, f, s)] = res
            fh.write(json.dumps({"param": k, "factor": f, "seed": s, "out": res}) + "\n")
            fh.flush()

    g0 = average([done[("__base__", 1.0, s)] for s in seeds])
    if "failed" in g0:
        raise SystemExit(f"base point failed: {g0['failed']}")
    S = {}
    for k in names:
        per_eps = []
        for e in EPS:
            gp = average([done[(k, 1 + e, s)] for s in seeds])
            gm = average([done[(k, 1 - e, s)] for s in seeds])
            if "failed" in gp or "failed" in gm:
                per_eps.append({j: np.nan for j in OUTPUTS})
                continue
            per_eps.append({j: (gp[j] - gm[j]) / (2 * e * g0[j]) if g0[j] else np.nan for j in OUTPUTS})
        S[k] = {j: float(np.nanmean([p[j] for p in per_eps])) for j in OUTPUTS}

    top = {j: sorted(names, key=lambda k: -abs(S[k][j]) if np.isfinite(S[k][j]) else 0)[:args.top] for j in OUTPUTS}
    freq = {k: sum(k in top[j] for j in OUTPUTS) for k in names}
    subset = sorted([k for k in names if freq[k] >= args.min_outputs], key=lambda k: -freq[k])

    json.dump({"base_outputs": g0, "sensitivities": S, "top_per_output": top, "frequency": freq,
               "selected": subset, "rule": f"top {args.top} in >= {args.min_outputs} of {len(OUTPUTS)} outputs"},
              open(out / "sa_summary.json", "w"), indent=2)
    json.dump(subset, open(out / "calibration_subset.json", "w"), indent=2)
    print("base outputs:", {k: round(v, 3) for k, v in g0.items() if isinstance(v, float)})
    print(f"selected ({len(subset)}):", ", ".join(f"{k}[{freq[k]}]" for k in subset))


if __name__ == "__main__":
    main()
