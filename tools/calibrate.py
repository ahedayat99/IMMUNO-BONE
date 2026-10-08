"""
Two-stage calibration (Methods 5.10.2): Latin hypercube screening, then differential evolution
refinement, for one of three nested search ranges around the base point.

  small  : base/3  .. base*3
  medium : base/10 .. base*10
  large  : params/calibration_bounds.json where available, otherwise base/100 .. base*100
Parameters taken from the literature are never moved outside their literature range
(LITERATURE_RANGES; single published values: +/-50 %, the rule used for COMMBINI's GA initialisation).
Sampling is log-uniform. Objective: macrophage loss (sum of squared relative errors of the four
immunofluorescence ratios), averaged over --seeds. Cytokines are not used for fitting.

DE (Methods Eq.): v = clip(x_best + F (x_a - x_b)), binomial crossover (CR), elitist selection.

Usage:
    python tools/calibrate.py --base params/option_B/sa_base.json --subset out/option_B/sa/calibration_subset.json \
        --range medium --out out/option_B/calib_medium --workers 32
"""
import argparse
import json
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

from evaluate import average, loss_macrophage, run_one, cytokine_directions

REPO = Path(__file__).resolve().parents[1]

# Published ranges (h^-1 or model units) for parameters that are fixed from the literature in Table S1.
LITERATURE_RANGES = {
    "k21": (2.08e-4, 2.08e-3),       # Trejo 2019 0.005-0.05/day
    "d1": (5.04e-3, 8.33e-3),        # Trejo 2019 0.121-0.2/day
    "d2": (6.79e-3, 8.33e-3),        # Trejo 2019 0.163-0.2/day
    "d_c2": (0.104, 0.193),          # Trejo 2019 2.5-4.632/day
    "d_c3": (0.417, 0.696),          # Trejo-GF 10/day; COMMBINI 1.16e-2/min
    "d_c4": (0.625, 1.875),          # Zhang 2021 / Geris 2008 30/day, +/-50 %
    "k_pm": (0.0208, 0.0421),        # Trejo 2020 0.5-1.01/day
    "a_pm": (0.01, 3.162),           # Trejo 2020 range
    "a_pm1": (5.0, 13.0),            # Trejo 2020 range
    "a_mb1": (0.1, 10.0),            # Trejo 2020 range (Kojouharov 2017)
    "k_e0": (0.125, 2.0), "k_e1": (0.125, 2.0), "k_e2": (0.125, 2.0),   # Trejo 2019 3-48/day
}
SINGLE_VALUE_FRACTION = 0.5          # +/-50 % for literature parameters with a single published value

# Optional per-job additions to LITERATURE_RANGES, e.g.
#   IMMUNO_LIT_EXTRA="d_c1:0.533:2.29;d_m_p:0.0144:0.0578"   (";" or "," separated)
# (used for round-3 variants; absent = unchanged behaviour).
import os as _os
for _item in [x for x in _os.environ.get("IMMUNO_LIT_EXTRA", "").replace(";", ",").split(",") if x]:
    _k, _lo, _hi = _item.split(":")
    LITERATURE_RANGES[_k] = (float(_lo), float(_hi))


def bounds_for(names, base, rng_name, original_calibrated, cal_bounds):
    lo, hi = [], []
    for k in names:
        b = base[k]
        if k in LITERATURE_RANGES:
            l, h = LITERATURE_RANGES[k]
        elif k not in original_calibrated:
            l, h = b * (1 - SINGLE_VALUE_FRACTION), b * (1 + SINGLE_VALUE_FRACTION)
        elif rng_name == "small":
            l, h = b / 3, b * 3
        elif rng_name == "medium":
            l, h = b / 10, b * 10
        else:
            l, h = (cal_bounds[k]["low"], cal_bounds[k]["high"]) if k in cal_bounds else (b / 100, b * 100)
        if k in original_calibrated and k in cal_bounds and rng_name != "large":
            l, h = max(l, cal_bounds[k]["low"]), min(h, cal_bounds[k]["high"])
        lo.append(l)
        hi.append(h)
    return np.log10(lo), np.log10(hi)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True)
    ap.add_argument("--subset", required=True)
    ap.add_argument("--range", choices=["small", "medium", "large"], required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--n_lhs", type=int, default=512)
    ap.add_argument("--lhs_seed_sim", type=int, default=1000)
    ap.add_argument("--seeds", default="1000,2000,3000", help="simulation seeds averaged in DE")
    ap.add_argument("--pop", type=int, default=64)
    ap.add_argument("--F", type=float, default=0.8)
    ap.add_argument("--CR", type=float, default=0.9)
    ap.add_argument("--max_gen", type=int, default=150)
    ap.add_argument("--patience", type=int, default=25, help="stop after this many generations without 0.5 %% improvement")
    ap.add_argument("--rng", type=int, default=7)
    args = ap.parse_args()

    base = json.load(open(args.base))
    names = json.load(open(args.subset))
    cal_bounds = {k: v for k, v in json.load(open(REPO / "params" / "calibration_bounds.json")).items() if not k.startswith("_")}
    original_calibrated = set(cal_bounds)
    lo, hi = bounds_for(names, base, args.range, original_calibrated, cal_bounds)
    rng = np.random.default_rng(args.rng)
    seeds = [int(s) for s in args.seeds.split(",")]
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    json.dump({"names": names, "log10_low": lo.tolist(), "log10_high": hi.tolist(), "range": args.range,
               "base": args.base}, open(out / "setup.json", "w"), indent=2)
    pool = ThreadPoolExecutor(max_workers=args.workers)

    def evaluate_many(X, sim_seeds):
        jobs = [(i, s) for i in range(len(X)) for s in sim_seeds]
        res = list(pool.map(lambda j: run_one(args.base, dict(zip(names, 10 ** X[j[0]])), j[1]), jobs))
        outs = [average(res[i * len(sim_seeds):(i + 1) * len(sim_seeds)]) for i in range(len(X))]
        return np.array([loss_macrophage(o) for o in outs]), outs

    # ---- Stage 1: Latin hypercube (one seed), resumable
    lhs_file = out / "lhs.jsonl"
    if lhs_file.exists():
        L = [json.loads(l) for l in open(lhs_file)]
    else:
        d = len(names)
        u = (rng.permuted(np.tile(np.arange(args.n_lhs), (d, 1)), axis=1).T + rng.random((args.n_lhs, d))) / args.n_lhs
        X = lo + u * (hi - lo)
        t0 = time.time()
        f, outs = evaluate_many(X, [args.lhs_seed_sim])
        L = [{"x": X[i].tolist(), "loss": float(f[i]), "out": outs[i]} for i in range(len(X))]
        with open(lhs_file, "w") as fh:
            for r in L:
                fh.write(json.dumps(r) + "\n")
        print(f"LHS done: {len(L)} samples in {time.time() - t0:.0f} s; best loss {min(r['loss'] for r in L):.4f}", flush=True)

    # ---- Stage 2: differential evolution (best/1/bin), initialised from the best LHS samples, resumable
    ck = out / "de_checkpoint.json"
    if ck.exists():
        c = json.load(open(ck))
        P, fP, gen, stale, hist = np.array(c["P"]), np.array(c["fP"]), c["gen"], c["stale"], c["history"]
    else:
        top = sorted(L, key=lambda r: r["loss"])[:args.pop]
        P = np.array([r["x"] for r in top])
        fP, _ = evaluate_many(P, seeds)          # re-score the initial population with all seeds
        gen, stale, hist = 0, 0, []
    best_prev = fP.min()
    while gen < args.max_gen and stale < args.patience:
        t0 = time.time()
        ib = int(np.argmin(fP))
        V = np.empty_like(P)
        for i in range(len(P)):
            a, b = rng.choice([j for j in range(len(P)) if j != i], 2, replace=False)
            v = np.clip(P[ib] + args.F * (P[a] - P[b]), lo, hi)
            cross = rng.random(len(names)) < args.CR
            cross[rng.integers(len(names))] = True
            V[i] = np.where(cross, v, P[i])
        fV, _ = evaluate_many(V, seeds)
        improved = fV < fP
        P[improved], fP[improved] = V[improved], fV[improved]
        gen += 1
        best = fP.min()
        stale = 0 if best < best_prev * 0.995 else stale + 1
        best_prev = min(best_prev, best)
        hist.append({"gen": gen, "best": float(best), "mean": float(np.mean(fP[fP < 1e3])), "sec": round(time.time() - t0, 1)})
        json.dump({"P": P.tolist(), "fP": fP.tolist(), "gen": gen, "stale": stale, "history": hist}, open(ck, "w"))
        print(f"gen {gen}: best {best:.4f} mean {hist[-1]['mean']:.4f} ({hist[-1]['sec']} s)", flush=True)

    # ---- Final: best parameter set, re-evaluated on 10 seeds (calibration + cytokine validation)
    ib = int(np.argmin(fP))
    best_params = dict(base)
    best_params.update(dict(zip(names, (10 ** P[ib]).tolist())))
    json.dump(best_params, open(out / "best_params.json", "w"), indent=4)
    final = [run_one(out / "best_params.json", {}, s) for s in range(1000, 11000, 1000)]
    avg = average(final)
    json.dump({"best_loss_de": float(fP[ib]), "final_10_seeds": avg, "loss_10_seeds": loss_macrophage(avg),
               "cytokines_rising": cytokine_directions(avg), "per_seed": final},
              open(out / "final_summary.json", "w"), indent=2)
    print("FINAL", json.dumps({k: (round(v, 3) if isinstance(v, float) else v) for k, v in avg.items()}), flush=True)


if __name__ == "__main__":
    main()
