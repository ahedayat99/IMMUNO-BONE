"""
Model-of-record analysis suite: baseline time courses + perturbation scenarios, 10 held-out seeds each
(same seeds in every scenario). Resumable: existing outputs are skipped.

Usage: python tools/run_final_suite.py --params params/round2/A2_P2__polish_05.json --out out/final --workers 8
"""
import argparse
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SEEDS = list(range(11000, 21000, 1000))


def scenarios(p):
    """Perturbation scenarios (Supplementary Table S2), as parameter overrides of the model of record."""
    msc0 = p.get("init_MSC_agents", 2)
    pmn0 = p.get("init_PMN_agents", 100)
    s = {
        "baseline": {},
        "pmn_half": {"init_PMN_agents": round(pmn0 * 0.5)},
        "pmn_10pct": {"init_PMN_agents": round(pmn0 * 0.1)},
        "pmn_double": {"init_PMN_agents": round(pmn0 * 2)},
        "krp_half": {"k_rp": p["k_rp"] * 0.5},
        "krp_double": {"k_rp": p["k_rp"] * 2},
        "m1_bias": {"k01": p["k01"] * 2, "k02": p["k02"] * 0.5},
        "m2_bias": {"k01": p["k01"] * 0.5, "k02": p["k02"] * 2},
        "msc_half": {"init_MSC_agents": max(1, round(msc0 * 0.5))},
        "msc_double": {"init_MSC_agents": round(msc0 * 2)},
        "kcm_half": {"k_cm": p["k_cm"] * 0.5},      # MSC VEGF secretion
        "kcm_double": {"k_cm": p["k_cm"] * 2},
        "kpm_half": {"k_pm": p["k_pm"] * 0.5},      # MSC proliferation rate (added)
        "kpm_double": {"k_pm": p["k_pm"] * 2},
        "klm_half": {"K_lm": p["K_lm"] * 0.5},
        "klm_double": {"K_lm": p["K_lm"] * 2},
        "debris_half": {"debris_coeff": 0.5},
        "debris_double": {"debris_coeff": 2.0},
        "debris_quad": {"debris_coeff": 4.0},
    }
    return s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--params", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--only", default="", help="comma-separated scenario names")
    args = ap.parse_args()
    p = json.load(open(args.params))
    out = Path(args.out)
    (out / "timecourse").mkdir(parents=True, exist_ok=True)
    S = scenarios(p)
    if args.only:
        S = {k: v for k, v in S.items() if k in args.only.split(",")}
    json.dump({k: v for k, v in S.items()}, open(out / "scenarios.json", "w"), indent=2)

    jobs = []
    for name, ov in S.items():
        for seed in SEEDS:
            f = out / "timecourse" / f"{name}_{seed}.json"
            if f.exists() and f.stat().st_size > 0:
                continue
            cmd = [sys.executable, str(REPO / "tools" / "timecourse.py"), "--params", args.params, "--seed", str(seed)]
            if name == "baseline":
                cmd += ["--snap", "0,24,48,72,96,120"]
            for k, v in ov.items():
                cmd += ["--set", f"{k}={v}"]
            jobs.append((f, cmd))
    print(f"{len(jobs)} runs to do", flush=True)

    def run(job):
        f, cmd = job
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            return f"FAILED {f.name}: {r.stderr[-300:]}"
        f.write_text(r.stdout)
        return f"ok {f.name}"

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for i, msg in enumerate(pool.map(run, jobs), 1):
            if msg.startswith("FAILED") or i % 20 == 0:
                print(f"[{i}/{len(jobs)}] {msg}", flush=True)
    print("done", flush=True)


if __name__ == "__main__":
    main()
