# Round 3 on SCC: TNF-α decay and neutrophil lifespan within literature ranges

Purpose: test whether the model of record can keep its fit when two biologically questionable values are held within literature ranges:
- **TNF-α decay `d_c1`**, currently at its search bound of 0.05 h⁻¹. Literature 0.533–2.29 h⁻¹ (Trejo et al. 2019, 12.79–55 day⁻¹). This addresses Reviewer comment 10.
- **Neutrophil apoptosis `d_m_p`**, currently 0.003 h⁻¹ (tissue half-life ≈ 10 days). New range 0.0144–0.099 h⁻¹ (half-life 7–48 h). The fast end is the neutrophil half-life of about 7 h (Summers et al. 2010, Trends Immunol 31:318). The slow end allows for delayed apoptosis in inflamed tissue; its citation is still to be confirmed by the authors. This targets the late neutrophil rise.

No model equations change. Only the bounds of these two parameters change.

## Install
Unzip `scc_round3_update.zip` in the repository root. It overwrites or adds:
- `tools/calibrate.py` — optional env var `IMMUNO_LIT_EXTRA` that appends literature ranges ("name:lo:hi;name:lo:hi"). Absent = unchanged.
- `tools/timecourse.py` — hourly totals, zones, parcel composition, capacity (analysis only).
- `scripts/*.py`, `tools/harness.py`, `tools/evaluate.py` — the exact local code used for the model of record (adds an `n_split` option, default 1 = identical results, and docstrings).
- `params/model_of_record.json` — the current model of record (layout 2, literature values clamped), used ONLY as a warm start.
- `params/selected_v2/loo_base_P2.json` — the base point: literature values + layout 2. The ±50 % bands of single literature values are therefore measured from the literature values, not from a previous fit.

**Verify:** `python tools/harness.py --seed 11000 --params params/model_of_record.json` must give
M1_M0_72 = 1.6454, M2_M0_72 = 0.1323, M1_M0_120 = 1.6987, M2_M0_120 = 1.4520. If not, STOP and report.

## Variants (existing `joint/` pipeline; run all three in parallel)
Common arguments, as for Stage A of round 2:
`--base params/selected_v2/loo_base_P2.json --tier wide --fit_cyt c2+c3+c4 --extra a02:0.0025:0.0075+a01:0.005:0.015 --warm_params params/model_of_record.json --warm out/joint/A2_P2+out/joint/A2_P2/polish`,
plus the same DE settings and 6-seed polishing as round 2. Warm-start points are clipped to the new bounds by the pipeline.

| Variant | Environment (pass with `qsub -v`) |
|---|---|
| R3_both | `IMMUNO_LIT_EXTRA="d_c1:0.533:2.29;d_m_p:0.0144:0.099"` |
| R3_tnf  | `IMMUNO_LIT_EXTRA="d_c1:0.533:2.29"` |
| R3_pmn  | `IMMUNO_LIT_EXTRA="d_m_p:0.0144:0.099"` |

Each coordinator **and** each polishing job must receive its variant's `IMMUNO_LIT_EXTRA`. Confirm it in each variant's `setup.json` (`log10_low/high` of `d_c1`, `d_m_p`).

## Selection and checks
1. Held-out selection on seeds 11000–20000 within each variant, using the same global objective as Stage A (c2+c3+c4).
2. For each variant's selected candidate, run on the 10 held-out seeds:
   - `tools/harness.py` (via the existing selection): macrophage ratios, MAPE, cytokine amount folds, perio_ctr_72 and perio_ctr_120;
   - `tools/timecourse.py`: report the mean PMN total at t = 0, 24, 48, 72, 96, 120 h, and total macrophages (M0+M1+M2) at 72 and 120 h.
3. Report the same items for `params/model_of_record.json` as the reference row.
4. For each selected candidate, report the final `d_c1` and `d_m_p`, and list every fitted parameter within 2 % (log scale) of a bound.

## Bring back
`round3_results.zip` containing:
- `out/joint/ROUND3_REPORT.md`;
- the three variant folders (setup.json, de_checkpoint.json, polish folders);
- the selected candidate JSONs with per-seed held-out results;
- the timecourse outputs;
- `logs/` and `RUN_NOTES.md`.
