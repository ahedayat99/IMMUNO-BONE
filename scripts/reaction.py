"""
Element-local reaction step (Methods 5.3-5.4; equation register E1-E13).

For one finite element e the coupled state is

    y = [D, c1, c2, c3, c4,  N_1, M0_1, M1_1, M2_1, C_1,  ...,  N_n, M0_n, M1_n, M2_n, C_n]

i.e. the element's debris and cytokine concentrations followed by the cell amounts of its n
resident parcels (each parcel is a mixed community of PMN, M0, M1, M2 and MSC).

Units: rates in h^-1. Cell amounts are cell-equivalents; capacities (PMN_max, M_max, K_lm)
and source terms are expressed per mean-sized element and scaled by a_e = A_e / mean(A),
so results do not depend on element size. Cytokines are in arbitrary concentration units.
"""
import numpy as np
from scipy.integrate import solve_ivp

from bone_healing_model_optimized import hypoxia_regulation
from parameters import DT_HOURS

N_FIELDS = 5   # D, c1, c2, c3, c4
N_CELLS = 5    # PMN, M0, M1, M2, MSC per parcel

# Debris engulfment weights applied to k_e (v1.0.1 update_debris_field; Table S1 kappa)
ENGULF_WEIGHT = {"PMN": 10.0, "M0": 5.0, "M1": 3.0, "M2": 2.0}


def element_rhs(t, y, p, n, a_e, EC, h):
    D, c1, c2, c3, c4 = np.maximum(y[:N_FIELDS], 0.0)
    cells = np.maximum(y[N_FIELDS:].reshape(n, N_CELLS), 0.0)
    N_p, M0_p, M1_p, M2_p, C_p = cells.T
    N, M0, M1, M2, C = cells.sum(axis=0)

    # E3: polarization transfer rates (Trejo et al. 2019; COMMBINI Eqs. 2-5)
    P01 = p["k01"] * c1 / (p["a01"] + c1)
    P02 = p["k02"] * c2 / (p["a02"] + c2)
    P12 = p["k12"] * c2 / (p["a_M1_to_M2"] + c2)
    P21 = p["k21"] * c1 / (p["a_M2_to_M1"] + c1)

    # E1, E2: element-level recruitment with density-based capacity, shared by the n parcels
    R_N = p["k_rp"] * D * a_e * max(0.0, 1.0 - N / (p["PMN_max"] * a_e))
    R_M = p["k_max"] * D * a_e * max(0.0, 1.0 - (M0 + M1 + M2) / (p["M_max"] * a_e))

    # E6, E7: MSC proliferation and differentiation
    A_m = p["k_pm"] * (p["a_pm"] ** 2 + p["a_pm1"] * c1) / (p["a_pm"] ** 2 + c1 ** 2)
    F1 = p["d_m"] * (p["a_mb1"] / (p["a_mb1"] + c1)) * (c3 / (p["a_mb"] + c3))

    dN = R_N / n - p["d_m_p"] * N_p
    dM0 = R_M / n - (P01 + P02 + p["d0"]) * M0_p
    dM1 = P01 * M0_p - (P12 + p["d1"]) * M1_p + P21 * M2_p
    dM2 = P02 * M0_p + P12 * M1_p - (P21 + p["d2"]) * M2_p
    dC = A_m * C_p * (1.0 - C / (p["K_lm"] * a_e)) - F1 * C_p   # E5 (-F1)

    # E8-E11: cytokine fields (sources per mean-sized element; decay applied only here)
    H1 = p["a12"] / (p["a12"] + c2 + c3)
    H2 = p["a22"] / (p["a22"] + c2)
    dc1 = H1 * (p["k0"] * D + (p["k1"] * M1 + p["k6"] * M0 + p["k7"] * N) / a_e) - p["d_c1"] * c1
    dc2 = H2 * (p["k8"] * M0 + p["k2"] * M2 + p["k3"] * C) / a_e - p["d_c2"] * c2
    dc3 = (p["k9"] * (M0 + M1) + p["k5"] * M2) / a_e - p["d_c3"] * c3
    dc4 = (h * (p["k_m0"] * M0 + p["k_m1"] * M1 + p["k_m2"] * M2 + p["k_cm"] * C) / a_e
           - p["d_c4"] * c4 - p["d_c4_ec"] * EC * c4)

    # E13: phagocytic debris clearance (Hill type II)
    engulf = (ENGULF_WEIGHT["PMN"] * p["k_e_pmn"] * N + ENGULF_WEIGHT["M0"] * p["k_e0"] * M0
              + ENGULF_WEIGHT["M1"] * p["k_e1"] * M1 + ENGULF_WEIGHT["M2"] * p["k_e2"] * M2)
    dD = -D / (p["a_ed"] + D) * engulf / a_e

    cell_rates = np.column_stack([dN, dM0, dM1, dM2, dC]).ravel()
    return np.concatenate([[dD, dc1, dc2, dc3, dc4], cell_rates])


def production_rates(fields, cells, p, a_e, h):
    """
    Cytokine secretion (production) rates of one element, as amounts per hour in mean-element units
    (the source terms of E8-E11 multiplied by a_e). Used to compare with dPCR transcript abundance,
    which reflects the rate of synthesis rather than accumulated protein.
    """
    D, c1, c2, c3, c4 = np.maximum(fields, 0.0)
    N, M0, M1, M2, C = np.maximum(cells, 0.0).sum(axis=0)
    H1 = p["a12"] / (p["a12"] + c2 + c3)
    H2 = p["a22"] / (p["a22"] + c2)
    return np.array([
        H1 * (p["k0"] * D * a_e + p["k1"] * M1 + p["k6"] * M0 + p["k7"] * N),
        H2 * (p["k8"] * M0 + p["k2"] * M2 + p["k3"] * C),
        p["k9"] * (M0 + M1) + p["k5"] * M2,
        h * (p["k_m0"] * M0 + p["k_m1"] * M1 + p["k_m2"] * M2 + p["k_cm"] * C),
    ])


def react_element(fields, cells, p, a_e, EC=0, PO2=0.0, rtol=1e-4, atol=1e-6, dt=DT_HOURS):
    """
    Integrate one element over dt hours (default DT_HOURS).
    fields: array [D, c1, c2, c3, c4]; cells: (n, 5) array of parcel amounts.
    Returns (new_fields, new_cells, n_clipped) with negatives clipped to zero.
    """
    n = cells.shape[0]
    h = hypoxia_regulation(PO2)
    y0 = np.concatenate([fields, cells.ravel()])
    sol = solve_ivp(element_rhs, (0.0, dt), y0, method="RK45", rtol=rtol, atol=atol,
                    args=(p, n, a_e, EC, h))
    y = sol.y[:, -1]
    n_clipped = int(np.sum(y < -atol))
    y = np.maximum(y, 0.0)
    return y[:N_FIELDS], y[N_FIELDS:].reshape(n, N_CELLS), n_clipped


def decay_unoccupied(c, p, dt=DT_HOURS):
    """Cytokines in elements without parcels only decay (exact solution over dt)."""
    rates = np.array([p["d_c1"], p["d_c2"], p["d_c3"], p["d_c4"]])
    return c * np.exp(-rates * dt)
