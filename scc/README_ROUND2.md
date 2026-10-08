# Round 2 on SCC: periosteal progenitor layouts + leave-one-cytokine-out cross-validation

Unzip `scc_round2_update.zip` in the repository root on SCC (`/projectnb/bonesgrp/ABM/IMMUNO-BONE_v2`).
It overwrites or adds only:
- `scripts/domain_model.py` — new `init_MSC_layout` switch. 0 = original (default, identical results),
  1 = progenitors in the periosteal cambium layer, 2 = half periosteal, half gap; total amount unchanged.
- `scripts/simple_domain_modular_v1.py` — accepts the two new optional keys.
- `tools/zone_explore.py` — zone-definition sensitivity (analysis only).
- `params/selected_v2/*.json` — base points:
  - `base_P1.json`, `base_P2.json`: current clamped selected set + layout 1 / 2 (Stage A).
  - `loo_base_P0/P1/P2.json`: pre-calibration Option-B base + layout (Stage B; never saw the cytokine targets).

Verify first: `python tools/harness.py --seed 11000 --params params/selected_v2/params_litclamped.json` must give
M1_M0_72 = 1.583, M2_M0_72 = 0.108, M1_M0_120 = 1.817, M2_M0_120 = 1.132 (single seed; layout 0 unchanged).

Common calibration settings (existing `joint/` pipeline, unchanged):
`--tier wide --extra a02:0.0025:0.0075+a01:0.005:0.015`. Under these settings:
- literature single values stay within ±50 %;
- published ranges are kept;
- the 13 originally calibrated parameters use their full bounds.

Use the same DE settings as the selected run, including 6-seed polishing of the best variant. Outputs go in **new**
variant folders. Do not mix them with the round-1 variants when selecting.

## Stage A — progenitor layout (run A1 and A2 in parallel)
| Variant | base | fit_cyt | warm starts |
|---|---|---|---|
| A1_P1 | params/selected_v2/base_P1.json | c2+c3+c4 | `--warm_params params/selected_v2/params_litclamped.json` + round-1 `p_x`, `x_c234` populations |
| A2_P2 | params/selected_v2/base_P2.json | c2+c3+c4 | same |

**Selection.** For each variant, re-run its 8 best distinct members on held-out seeds 11000–20000 (as `joint/select_final.py`,
restricted to the A1/A2 variants). Report the best candidate of each, and of the current clamped set (layout 0):
- macrophage MAPE (targets 1.60 / 0.14 / 1.69 / 1.37);
- the four cytokine amount folds;
- `perio_ctr_72` and `perio_ctr_120`.

**Pre-registered layout choice. Use the global objective ONLY, never the zone values:**
- Lowest held-out c2+c3+c4 objective wins.
- If P1 and P2 are within 10 % of each other, choose P2.
- If neither beats layout 0 by more than 10 %, layout 0 stays and the zone result is reported as is.

## Stage B — leave-one-cytokine-out cross-validation (on the chosen layout L; run B1–B3 in parallel)
| Variant | base | fit_cyt | omitted (predicted) |
|---|---|---|---|
| B1_noIL10 | params/selected_v2/loo_base_L.json | c3+c4 | IL-10 |
| B2_noTGF | params/selected_v2/loo_base_L.json | c2+c4 | TGF-β |
| B3_noVEGF | params/selected_v2/loo_base_L.json | c2+c3 | VEGFA |

Additional settings for every Stage B variant:
- `--w_dir 0`, so the omitted cytokine does not enter the loss in any form.
- Warm starts ONLY from the macrophage-only runs (`out/option_B/calib_small`, `out/option_B/calib_medium`).
  Do **not** warm-start from any round-1 or Stage A population, nor from params_litclamped.json.
  Those have seen all cytokine targets.

Select each variant's best candidate on held-out seeds 11000–20000, exactly as in Stage A.
Report, for each variant:
- the macrophage MAPE;
- the fitted cytokines' errors;
- **the predicted (omitted) cytokine's fold and error**.

## Bring back
- `out/joint/A1_P1`, `A2_P2`, `B1_noIL10`, `B2_noTGF`, `B3_noVEGF`: setup.json, de_checkpoint.json, and any polish folders.
- The held-out selection outputs: candidate JSONs and per-seed results of the selected candidates.
- `logs/`.
