"""
Bone healing model equations and related functions - OPTIMIZED VERSION
"""


# Move hypoxia_regulation outside to avoid redefining it every call
def hypoxia_regulation(oxygen):
    """Calculate hypoxia regulation factor based on oxygen level."""
    if oxygen < 1.0:
        return 6.0
    elif oxygen <= 20.0:
        return 1.0 + 5.0 * ((20.0 - oxygen) / 19.0) ** 3
    else:
        return 1.0


def full_bone_healing_model(t, variables, params, D, EC=0, PO2=0.0):
    """
    Define the bone healing model equations - OPTIMIZED.

    :param t: Time parameter
    :param variables: State variables [PMN, M0, M1, M2, c1, c2, c3, Cm, c4]
    :param params: Dictionary of model parameters
    :param D: Debris concentration
    :param EC: Endothelial cell count (default 0)
    :param PO2: Oxygen partial pressure (default 0.0)
    :return: List of derivatives for each state variable
    """
    PMN, M0, M1, M2, c1, c2, c3, Cm, c4 = variables

    # Extract parameters - still needs optimization via caching at agent level
    k_e0, k_e1, k_e2, k_e_pmn = params["k_e0"], params["k_e1"], params["k_e2"], params["k_e_pmn"]
    a_mb1, a_mb = params["a_mb1"], params["a_mb"]
    d0, d1, d2, d_m_p = params["d0"], params["d1"], params["d2"], params["d_m_p"]
    k0, k1, k2, k3, k5, k6, k7, k8, k9, d_c1, d_c2, d_c3 = params["k0"], params["k1"], params["k2"], params["k3"], params["k5"], params["k6"], params["k7"], params["k8"], params["k9"], params["d_c1"], params["d_c2"], params["d_c3"]
    k_pm, a_pm, a_pm1 = params["k_pm"], params["a_pm"], params["a_pm1"]
    a_ed, k_max, M_max, PMN_max = params["a_ed"], params["k_max"], params["M_max"], params["PMN_max"]
    k_M0_to_M1, a_M0_to_M1 = params["k01"], params["a01"]
    k_M0_to_M2, a_M0_to_M2 = params["k02"], params["a02"]
    k_M1_to_M2, a_M1_to_M2 = params["k12"], params["a_M1_to_M2"]
    k_M2_to_M1, a_M2_to_M1 = params["k21"], params["a_M2_to_M1"]
    k_rp, d_m = params["k_rp"], params["d_m"]
    k_m2, k_m1, k_m0, k_cm = params["k_m2"], params["k_m1"], params["k_m0"], params["k_cm"]
    d_c4, d_c4_ec = params["d_c4"], params["d_c4_ec"]
    a12, a22, a33 = params["a12"], params["a22"], params["a33"]
    K_lm = params["K_lm"]

    # Polarization probabilities
    P_M0_to_M1 = k_M0_to_M1 * (c1 / (a_M0_to_M1 + c1))
    P_M0_to_M2 = k_M0_to_M2 * (c2 / (a_M0_to_M2 + c2))
    P_M1_to_M2 = k_M1_to_M2 * (c2 / (a_M1_to_M2 + c2))
    P_M2_to_M1 = k_M2_to_M1 * (c1 / (a_M2_to_M1 + c1))

    # Recruitment and regulation
    R_D = D / (a_ed + D)
    R_PMN = k_rp * (1 - (PMN / PMN_max)) * D
    R_M = k_max * (1 - (M0 + M1 + M2) / M_max) * D
    H1 = a12 / (a12 + c2 + c3)
    H2 = a22 / (a22 + c2)
    H3 = c3 / (a33 + c3)
    A_m = k_pm * (a_pm**2 + a_pm1 * c1) / (a_pm**2 + c1**2)
    F1 = d_m * (a_mb1 / (a_mb1 + c1)) * (c3 / (a_mb + c3))

    # Derivatives
    dPMN_dt = R_PMN - d_m_p * PMN
    dM0_dt = R_M - P_M0_to_M1 * M0 - P_M0_to_M2 * M0 - d0 * M0
    dM1_dt = P_M0_to_M1 * M0 - d1 * M1 - P_M1_to_M2 * M1 + P_M2_to_M1 * M2
    dM2_dt = P_M0_to_M2 * M0 - d2 * M2 + P_M1_to_M2 * M1 - P_M2_to_M1 * M2
    dc1_dt = H1 * (k0 * D + k1 * M1 + k6 * M0 + k7 * PMN) - d_c1 * c1
    dc2_dt = H2 * (k8 * M0 + k2 * M2 + k3 * Cm) - d_c2 * c2
    dc3_dt = (k9 * (M0 + M1) + k5 * M2) - d_c3 * c3

    # VEGF dynamics with hypoxia regulation
    hypoxia_factor = hypoxia_regulation(PO2)
    dc4_dt = hypoxia_factor * (
        k_m2 * M2 +
        k_m1 * M1 +
        k_m0 * M0 +
        k_cm * Cm
    ) - d_c4 * c4 - d_c4_ec * EC
    dCm_dt = A_m * Cm * (1 - Cm / K_lm) + F1 * Cm

    return [dPMN_dt, dM0_dt, dM1_dt, dM2_dt, dc1_dt, dc2_dt, dc3_dt, dCm_dt, dc4_dt]
