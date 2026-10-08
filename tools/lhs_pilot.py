"""
Latin-hypercube screen of the 13 calibrated parameters (log-uniform within
params/calibration_bounds.json); all other parameters fixed at the chosen option.

Each sample is run with tools/harness.py (one seed) and appended to a JSONL file together
with two losses (sum of squared relative errors):
  loss_M   : four immunofluorescence macrophage ratios only
  loss_all : the four ratios plus the four dPCR cytokine fold-changes

Usage:
    python tools/lhs_pilot.py --option B --n 400 --seed 1000 --workers 8 --out out/option_B/pilot_lhs.jsonl
"""
import argparse
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
TARGETS_M = {"M1_M0_72": 1.620, "M2_M0_72": 0.145, "M1_M0_120": 1.667, "M2_M0_120": 1.367}
TARGETS_C = {"c1": 12.85, "c2": 6.20, "c3": 8.97, "c4": 5.54}


def lhs_log(bounds, n, rng):
    names = list(bounds)
    d = len(names)
    u = (rng.permuted(np.tile(np.arange(n), (d, 1)), axis=1).T + rng.random((n, d))) / n
    lo = np.log10([bounds[k]["low"] for k in names])
    hi = np.log10([bounds[k]["high"] for k in names])
    return names, 10 ** (lo + u * (hi - lo))


def losses(rec):
    def rel(p, t):
        return ((p - t) / t) ** 2 if p is not None and np.isfinite(p) else 1e3
    lm = sum(rel(rec["calibration"][k], t) for k, t in TARGETS_M.items())
    lc = sum(rel(rec["cytokine_fold"][k], t) for k, t in TARGETS_C.items())
    return lm, lm + lc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--option", choices=["A", "B"], required=True)
    ap.add_argument("--n", type=int, default=400)
    ap.add_argument("--seed", type=int, default=1000, help="simulation seed (same for all samples)")
    ap.add_argument("--lhs_seed", type=int, default=42)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--out", required=True)
    ap.add_argument("--max_parcels", type=int, default=3000)
    args = ap.parse_args()

    bounds = {k: v for k, v in json.load(open(REPO / "params" / "calibration_bounds.json")).items()
              if not k.startswith("_")}
    names, X = lhs_log(bounds, args.n, np.random.default_rng(args.lhs_seed))
    params_file = REPO / "params" / f"option_{args.option}" / "params.json"
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    def run(i):
        cmd = [sys.executable, str(REPO / "tools" / "harness.py"), "--seed", str(args.seed),
               "--params", str(params_file), "--label", f"lhs{i}", "--max_parcels", str(args.max_parcels)]
        for k, v in zip(names, X[i]):
            cmd += ["--set", f"{k}={v:.8g}"]
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
        if r.returncode != 0:
            return {"i": i, "params": dict(zip(names, X[i].tolist())), "error": r.stderr[-500:]}
        rec = json.loads(r.stdout.strip().splitlines()[-1])
        if rec.get("aborted"):
            return {"i": i, "params": dict(zip(names, X[i].tolist())), "aborted": True,
                    "hour": rec["hour"], "n_parcels": rec["n_parcels"], "loss_M": 1e3, "loss_all": 1e3}
        lm, la = losses(rec)
        return {"i": i, "params": dict(zip(names, X[i].tolist())), "loss_M": lm, "loss_all": la,
                "calibration": rec["calibration"], "cytokine_fold": rec["cytokine_fold"],
                "parcels_120": rec["snapshots"]["120"]["n_parcels"],
                "zones_120": rec["snapshots"]["120"]["zones"], "zones_72": rec["snapshots"]["72"]["zones"],
                "totals_72": rec["snapshots"]["72"]["totals"], "totals_120": rec["snapshots"]["120"]["totals"],
                "elapsed_s": rec["elapsed_s"]}

    done = set()
    if out.exists():
        done = {json.loads(l)["i"] for l in open(out) if l.strip()}
    todo = [i for i in range(args.n) if i not in done]
    print(f"{len(done)} samples already done, running {len(todo)}")
    with ThreadPoolExecutor(max_workers=args.workers) as pool, open(out, "a") as f:
        for fut in as_completed([pool.submit(run, i) for i in todo]):
            f.write(json.dumps(fut.result()) + "\n")
            f.flush()
    print(f"wrote {args.n} samples to {out}")


if __name__ == "__main__":
    main()
