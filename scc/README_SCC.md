# Running the v2 sensitivity analysis and calibration on BU SCC

## 1. Upload
Copy the whole `IMMUNO-BONE_v2/` folder (about 3 MB) to SCC, for example to `/projectnb/bonesgrp/ABM/IMMUNO-BONE_v2`.

## 2. One-time environment
```bash
cd /projectnb/bonesgrp/ABM/IMMUNO-BONE_v2
bash scc/setup_env.sh          # creates conda env "immuno_v2" with mesa==2.4.0 (Mesa 3 will NOT work)
```

## 3. Quick check on a login node (about 30 s)
```bash
module load miniconda && conda activate immuno_v2
python tools/harness.py --seed 1000 --params params/option_B/sa_base.json | python -c "import json,sys; d=json.load(sys.stdin); print(d['calibration'])"
# expected (same seed, same code): M1_M0_72 ~1.33, M2_M0_72 ~0.138, M1_M0_120 ~1.81, M2_M0_120 ~0.557
```

## 4. Submit (Option B)
```bash
bash scc/submit_pipeline.sh
```
This submits:
- **One sensitivity-analysis job**, `immuno_sa`: 957 runs, about 15–30 min on 32 cores.
- **Three calibration jobs**, `immuno_calib` for the small, medium and large ranges. They start automatically when the sensitivity job finishes. Each runs an LHS of 512 samples, then DE with population 64 and 3 seeds per evaluation, stopping early after 25 stale generations. Each takes a few hours.

Monitor with `qstat -u $USER`; logs are in `logs/`. All steps are resumable: re-submitting continues from `lhs.jsonl` / `de_checkpoint.json` / `raw_runs.jsonl`.

## 5. Outputs to copy back
```
out/option_B/sa/                 sa_summary.json, calibration_subset.json, raw_runs.jsonl
out/option_B/calib_small/        setup.json, lhs.jsonl, de_checkpoint.json, best_params.json, final_summary.json
out/option_B/calib_medium/       (same)
out/option_B/calib_large/        (same)
logs/
```

## 6. Fallback: Option A (only if Option B cannot reach the targets)
```bash
python tools/lhs_pilot.py --option A --n 400 --workers 31 --out out/option_A/pilot_lhs.jsonl
python tools/make_sa_base.py --option A --pilot out/option_A/pilot_lhs.jsonl
OPTION=A bash scc/submit_pipeline.sh
```
