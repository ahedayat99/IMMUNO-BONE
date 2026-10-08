"""
Parameter bookkeeping shared by the model and the tools.

Time base
---------
The ODE system is integrated over one global step of DT_HOURS = 1 h, so every rate
constant is applied at its stated value in h^-1, as listed in Supplementary Table S1.

Code versions up to v1.0.1 integrated the ODEs over only 0.1 time units per 1-h step,
i.e. they applied every ODE rate constant at one tenth of its stated value.
v1_to_per_hour() reproduces that v1.0.1 behaviour under the 1-h step. It is used only
for the regression test against the archived v1.0.1 outputs
(params/v1_behaviour_equivalent_per_hour.json), never for results.
"""

DT_HOURS = 1.0

# Parameters that multiply a right-hand-side term of full_bone_healing_model() (h^-1).
ODE_RATE_KEYS = (
    # polarization maximum rates
    "k01", "k02", "k12", "k21",
    # macrophage and PMN removal
    "d0", "d1", "d2", "d_m_p",
    # recruitment
    "k_rp", "k_max",
    # cytokine production
    "k0", "k1", "k2", "k3", "k5", "k6", "k7", "k8", "k9",
    "k_m0", "k_m1", "k_m2", "k_cm",
    # cytokine decay and endothelial VEGF uptake
    "d_c1", "d_c2", "d_c3", "d_c4", "d_c4_ec",
    # MSC proliferation and differentiation
    "k_pm", "d_m",
)

V1_ODE_TIME_UNIT_H = 0.1


def v1_to_per_hour(params_v1):
    """Convert a v1.0.x parameter set (ODE rates per 0.1 h) to per-hour rates."""
    converted = dict(params_v1)
    for key in ODE_RATE_KEYS:
        converted[key] = params_v1[key] * V1_ODE_TIME_UNIT_H
    return converted
