#!/usr/bin/env python3
"""
Main simulation file for bone healing model.
Runs the simulation and generates output profiles.

Usage:
  python scripts/simple_domain_modular_v1.py --node data/node_elements.txt
  python scripts/simple_domain_modular_v1.py --node data/node_elements.txt --verbose
  python scripts/simple_domain_modular_v1.py --node data/node_elements.txt --params_json test_run/best_params.json
  python scripts/simple_domain_modular_v1.py --node data/node_elements.txt --output_csv ./output/my_run.csv
"""

import time
import csv
import argparse
import sys
import json
from pathlib import Path

from utils import parse_nodes_elements
from domain_model import DomainModel


# =====================================================
# Simulation Functions
# =====================================================

def run_to_targets(nodes, elements, params, targets=(72, 120), disable_EC=True):
    model = DomainModel(nodes, elements, params)

    if hasattr(model, "debug"):
        model.debug = False
    if hasattr(model, "enable_EC") and disable_EC:
        model.enable_EC = False

    targets = sorted(set(targets))
    out = {}

    while model.time < 120:
        model.step()
        if model.time in targets:
            out[model.time] = model.snapshot_outputs()

    return out


def run_and_save_all_outputs(
    nodes,
    elements,
    params,
    max_hours=120,
    disable_EC=True,
    output_file="simulation_outputs.csv",
    verbose=False
):
    from element_agent_optimized import ElementAgent

    model = DomainModel(nodes, elements, params)

    if hasattr(model, "debug"):
        model.debug = False
    if hasattr(model, "enable_EC") and disable_EC:
        model.enable_EC = False

    all_outputs = []

    for hour in range(1, max_hours + 1):
        model.step()
        snapshot = model.snapshot_outputs()

        total_debris = sum(model.debris_field.values())
        num_agents = len([
            a for a in model.schedule.agents
            if isinstance(a, ElementAgent)
        ])

        output_row = {
            "hour": hour,
            "total_PMN": snapshot["counts"]["PMN"],
            "total_M0": snapshot["counts"]["M0"],
            "total_M1": snapshot["counts"]["M1"],
            "total_M2": snapshot["counts"]["M2"],
            "total_MSC": snapshot["counts"]["MSC"],
            "total_EC": snapshot["counts"]["EC"],
            "total_c1": snapshot["cytokines"]["c1"],
            "total_c2": snapshot["cytokines"]["c2"],
            "total_c3": snapshot["cytokines"]["c3"],
            "total_c4": snapshot["cytokines"]["c4"],
            "total_debris": total_debris,
            "num_agents": num_agents
        }

        all_outputs.append(output_row)

        if verbose:
            print(
                f"[Hour {hour:3d}] "
                f"PMN={output_row['total_PMN']:.2f} | "
                f"M0={output_row['total_M0']:.2f} | "
                f"M1={output_row['total_M1']:.2f} | "
                f"M2={output_row['total_M2']:.2f} | "
                f"MSC={output_row['total_MSC']:.2f} | "
                f"EC={output_row['total_EC']:.2f} | "
                f"c1={output_row['total_c1']:.2f} | "
                f"c2={output_row['total_c2']:.2f} | "
                f"c3={output_row['total_c3']:.2f} | "
                f"c4={output_row['total_c4']:.2f} | "
                f"Debris={output_row['total_debris']:.2f} | "
                f"Agents={output_row['num_agents']}"
            )

    output_path = Path(output_file).resolve()

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", newline="") as f:
        fieldnames = [
            "hour", "total_PMN", "total_M0", "total_M1", "total_M2",
            "total_MSC", "total_EC",
            "total_c1", "total_c2", "total_c3", "total_c4",
            "total_debris", "num_agents"
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_outputs)

    print(f"Saved all outputs to: {output_path}")

    return all_outputs


# =====================================================
# Default Parameters
# =====================================================

DEFAULT_PARAMS = {
    "k_e0": 3.33e-2 * 60,
    "k_e1": 3.33e-2 * 60,
    "k_e2": 3.33e-2 * 60,
    "k_e_pmn": 3.33e-3 * 60,
    "k12": 5.8e-4 * 600,
    "k21": 3.47e-6 * 60,
    "d0": 1.08e-4 * 1000,
    "d1": 1.38e-4 * 60,
    "d2": 1.38e-4 * 60,
    "d_m_p": 7.00e-5 * 600,
    "k0": 5.9e-6 * 600,
    "k1": 5.76e-6 * 600,
    "k2": 2.58e-6 * 6000,
    "k3": 5.55e-6 * 6000,
    "k5": 2.67e-7 * 6000,
    "k6": 5.78e-7 * 600,
    "k7": 4.61e-7 * 600,
    "k8": 1.27e-5 * 6000,
    "k9": 3.13e-5 * 6000,
    "d_c1": 8.9e-4 * 60,
    "d_c2": 3.22e-4 * 60,
    "d_c3": 1.16e-4 * 60,
    "k_pm": 3.47e-4 * 250,
    "a_pm": 3.162,
    "a_pm1": 13,
    "d_m": 6.94e-4 * 6000,
    "a_mb1": 5,
    "a_mb": 1,
    "a_ed": 4.71e2,
    "k_max": 0.0084,
    "M_max": 2.5,
    "PMN_max": 10,
    "k01": 0.42519417,
    "k02": 0.04915,
    "a01": 0.01,
    "a02": 0.061165,
    "a12": 0.025,
    "a22": 0.6,
    "a33": 1.13585,
    "K_lm": 15,
    "k_rp": 1.04e-5,
    "a_M1_to_M2": 0.1,
    "a_M2_to_M1": 0.5,
    "diff_c1": 1.8e-5 * 60,
    "diff_c2": 1.8e-5 * 60,
    "diff_c3": 1.56e-5 * 60,
    "diff_c4": 1.5e-5 * 60,
    "d_c4": 0.0005,
    "d_c4_ec": 0.05,
    "k_m2": 3.0e-2,
    "k_m1": 1.5e-3,
    "k_m0": 5.0e-4,
    "k_cm": 2.0e-3
}


# =====================================================
# Params loading + strict validation
# =====================================================

def load_params_or_die(params_json_path: str | None, defaults: dict) -> dict:
    """
    If params_json_path is None -> return defaults.
    Else -> load JSON, require keys match defaults exactly (no missing, no extras),
    except for OPTIONAL initial-state keys which use defaults if absent.
    Also enforces numeric values (int/float).
    """
    OPTIONAL_INIT_KEYS = {
        # Existing optional initial-state keys
        "init_M0": 0.0,
        "init_M1": 0.0,
        "init_M2": 0.0,
        "init_c1": 0.050157,
        "init_c2": 0.0,
        "init_c3": 0.0,
        "init_c4": 0.02,

        # NEW optional agent-count keys (so schema check won't reject them)
        "init_PMN_agents": 100,
        "init_M0_agents": 1,
        "init_MSC_agents": 2,

        # NEW debris scaling coefficient
        "debris_coeff": 1.0,
    }

    if params_json_path is None:
        merged = dict(defaults)
        for k, v in OPTIONAL_INIT_KEYS.items():
            merged.setdefault(k, v)
        return merged

    p = Path(params_json_path)
    if not p.exists():
        print(f"\nERROR: params JSON file not found: {p}\n", file=sys.stderr)
        sys.exit(1)

    try:
        with open(p, "r", encoding="utf-8") as f:
            loaded = json.load(f)
    except json.JSONDecodeError as e:
        print(f"\nERROR: Failed to parse JSON in {p}\n  {e}\n", file=sys.stderr)
        sys.exit(1)

    if not isinstance(loaded, dict):
        print(
            f"\nERROR: Params JSON must be a JSON object/dict at the top-level.\n"
            f"Got type: {type(loaded).__name__}\n",
            file=sys.stderr
        )
        sys.exit(1)

    default_keys = set(defaults.keys())
    optional_keys = set(OPTIONAL_INIT_KEYS.keys())
    loaded_keys = set(loaded.keys())

    # Required keys = default keys (must all be present)
    # Optional keys = init_* keys (fill in defaults if absent)
    # Extra keys = anything not in required or optional
    missing = sorted(default_keys - loaded_keys)
    extra = sorted(loaded_keys - default_keys - optional_keys)

    if missing or extra:
        print("\nERROR: Invalid params JSON schema.\n", file=sys.stderr)
        if missing:
            print("Missing required parameter(s):", file=sys.stderr)
            for k in missing:
                print(f"  - {k}", file=sys.stderr)
        if extra:
            print("\nUnexpected extra parameter(s):", file=sys.stderr)
            for k in extra:
                print(f"  - {k}", file=sys.stderr)

        print(
            "\nFix: Your JSON must contain exactly these required keys:\n"
            + ", ".join(sorted(default_keys))
            + "\n\nOptionally, you may also include these initial-state keys:\n"
            + ", ".join(sorted(optional_keys))
            + "\n",
            file=sys.stderr
        )
        sys.exit(1)

    # Type check values
    non_numeric = []
    for k, v in loaded.items():
        if not isinstance(v, (int, float)):
            non_numeric.append((k, type(v).__name__))

    if non_numeric:
        print("\nERROR: All parameter values must be numeric (int/float).\n", file=sys.stderr)
        print("Non-numeric values found:", file=sys.stderr)
        for k, tname in non_numeric:
            print(f"  - {k}: {tname}", file=sys.stderr)
        sys.exit(1)

    # Fill in optional keys with defaults if not provided
    for k, v in OPTIONAL_INIT_KEYS.items():
        loaded.setdefault(k, v)

    return loaded


# =====================================================
# Argument Parser
# =====================================================

def parse_args():
    parser = argparse.ArgumentParser(description="Run bone healing simulation")

    parser.add_argument(
        "--node",
        required=True,
        help="Path to node_elements.txt file"
    )

    parser.add_argument(
        "--params_json",
        default=None,
        help=(
            "Optional path to JSON file containing model parameters. "
            "If provided, it must contain EXACTLY the same keys as DEFAULT_PARAMS "
            "(no missing keys, no extra keys)."
        )
    )

    parser.add_argument(
        "--output_csv",
        default=None,
        help=(
            "Optional path to output CSV file. "
            "If not provided, will overwrite the default: simulation_outputs_1_to_120_opt.csv "
            "(in the current working directory)."
        )
    )

    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print simulation outputs at every hour"
    )

    return parser.parse_args()


# =====================================================
# Main
# =====================================================

def main():
    args = parse_args()

    node_file = Path(args.node)
    if not node_file.exists():
        print("\nERROR: node_elements.txt file not found.\n", file=sys.stderr)
        print("Please provide the input file using:\n", file=sys.stderr)
        print("  python main_simulation.py --node data/node_elements.txt\n", file=sys.stderr)
        sys.exit(1)

    params = load_params_or_die(args.params_json, DEFAULT_PARAMS)

    print(f"Using node file: {node_file}")
    if args.params_json is None:
        print("Using DEFAULT_PARAMS (built-in).")
    else:
        print(f"Using params from JSON: {Path(args.params_json).resolve()}")

    nodes, elements = parse_nodes_elements(node_file)

    output_csv_path = args.output_csv if args.output_csv is not None else "simulation_outputs_1_to_120_opt.csv"

    print("Running simulation from hour 1 to 120...")
    start_time = time.perf_counter()

    all_outputs = run_and_save_all_outputs(
        nodes,
        elements,
        params,
        max_hours=120,
        disable_EC=True,
        output_file=output_csv_path,
        verbose=args.verbose
    )

    total_time = time.perf_counter() - start_time

    print(f"\nSimulation completed in {total_time:.2f} seconds")
    print(f"Average time per hour: {total_time/120:.4f} seconds")

    print("\nFirst hour summary:")
    print(f"  PMN: {all_outputs[0]['total_PMN']:.2f}")
    print(f"  M0: {all_outputs[0]['total_M0']:.2f}")
    print(f"  M1: {all_outputs[0]['total_M1']:.2f}")
    print(f"  M2: {all_outputs[0]['total_M2']:.2f}")
    print(f"  MSC: {all_outputs[0]['total_MSC']:.2f}")
    print(f"  Debris: {all_outputs[0]['total_debris']:.2f}")

    print("\nLast hour summary:")
    print(f"  PMN: {all_outputs[-1]['total_PMN']:.2f}")
    print(f"  M0: {all_outputs[-1]['total_M0']:.2f}")
    print(f"  M1: {all_outputs[-1]['total_M1']:.2f}")
    print(f"  M2: {all_outputs[-1]['total_M2']:.2f}")
    print(f"  MSC: {all_outputs[-1]['total_MSC']:.2f}")
    print(f"  Debris: {all_outputs[-1]['total_debris']:.2f}")


if __name__ == "__main__":
    main()