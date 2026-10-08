# Option B — non-calibrated parameters at their literature (source) values

Identical to Option A except for the nine non-calibrated parameters below, which are set to the value of the cited source (converted to h^-1 where needed). See papers/final/reviewer_response/parameter_provenance_audit.md.
The 13 calibrated parameters start from the v1.0.1 values and are re-fitted.

| Parameter | Option A (Table S1) | Option B | B/A | Source |
|---|---|---|---|---|
| d_c2 | 0.01932 | 0.1932 | 10 | COMMBINI Suppl. Table 2: IL-10 decay 3.22e-3 min^-1 (x60); = Trejo 2019 upper bound 4.632/day |
| d_c3 | 0.00696 | 0.696 | 100 | COMMBINI Suppl. Table 2: TGF-b decay 1.16e-2 min^-1 (x60) |
| k_rp | 1.04e-05 | 0.000624 | 60 | Trejo 2019 k_max 0.015/day = 1.04e-5 min^-1 (x60); Table S1 omitted the x60 |
| d_m_p | 0.042 | 0.0042 | 0.1 | COMMBINI Suppl. Table 1: PMN apoptosis 7.00e-5 min^-1 (x60) |
| k_pm | 0.08675 | 0.02082 | 0.24 | Trejo 2019 / Trejo 2020: MSC proliferation 0.5/day = 3.47e-4 min^-1 (x60) |
| a02 | 0.061165 | 0.005 | 0.0817 | Trejo 2019: half-saturation of c2 to activate M2, 0.005 ng/mL |
| a22 | 0.6 | 0.1 | 0.167 | Trejo 2019: effectiveness of c2 inhibition of c2 synthesis, 0.1 ng/mL |
| a_mb | 1 | 10 | 10 | Trejo et al. 2019 (AIP, growth factors): effectiveness of c3 enhancing Cm differentiation, 10 ng/mL |
| d_c4 | 0.0005 | 1.25 | 2500 | Zhang et al. 2021 J Biomech, Suppl. Table A2 (from Geris et al. 2008): VEGF decay d_gv = 30 day^-1 |

Diffusion coefficients are set in mm^2 h^-1 when the finite-volume transport step is implemented (task T7): COMMBINI values for TNF-a, IL-10, TGF-b (0.108, 0.108, 0.0936) and Zhang et al. 2021 / Geris et al. 2008 for VEGF (D_gv = 6.125e-6 m^2 day^-1 = 0.255).
