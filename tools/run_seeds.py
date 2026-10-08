"""
Run tools/harness.py for several seeds in parallel (one process per run) and
write one JSON record per line.

Usage:
    python tools/run_seeds.py --out tests/golden/v1_baseline.jsonl --seeds 1000:10000:1000 \
        [--scripts scripts] [--params data/input_params.json] [--set K_lm=10.06] [--workers 5]
"""
import argparse
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def parse_seeds(spec):
    if ":" in spec:
        start, stop, step = (int(v) for v in spec.split(":"))
        return list(range(start, stop + 1, step))
    return [int(v) for v in spec.split(",")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--seeds", default="1000:10000:1000")
    ap.add_argument("--scripts", default=str(REPO / "scripts"))
    ap.add_argument("--params", default=str(REPO / "data" / "input_params.json"))
    ap.add_argument("--set", action="append", default=[])
    ap.add_argument("--label", default="")
    ap.add_argument("--workers", type=int, default=5)
    args = ap.parse_args()

    def run(seed):
        cmd = [sys.executable, str(REPO / "tools" / "harness.py"), "--seed", str(seed),
               "--scripts", args.scripts, "--params", args.params, "--label", args.label]
        for s in args.set:
            cmd += ["--set", s]
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            raise RuntimeError(f"seed {seed} failed:\n{res.stderr[-2000:]}")
        return res.stdout.strip().splitlines()[-1]

    seeds = parse_seeds(args.seeds)
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        lines = list(pool.map(run, seeds))
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text("\n".join(lines) + "\n")
    print(f"wrote {len(lines)} runs to {args.out}")


if __name__ == "__main__":
    main()
