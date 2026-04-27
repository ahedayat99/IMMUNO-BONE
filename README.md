# Bone Healing Agent-Based Model (ABM)

This is a repository for a code based IMMUNO-BONE model

An agent-based model of the early inflammatory and regenerative phases of bone fracture healing. The model simulates spatial interactions among key cell populations — neutrophils (PMN), macrophages (M0/M1/M2), mesenchymal stem cells (MSC), and endothelial cells (EC) — along with pro- and anti-inflammatory cytokines, across a finite-element mesh derived from a bone callus geometry.

---

## Repository Structure

```
.
├── data/
│   ├── node_elements.txt   # FE mesh geometry (nodes + elements) of the bone callus
│   └── input_params.json   # Calibrated model parameters
├── scripts/
│   ├── simple_domain_modular_v1.py      # Main entry point — run this
│   ├── domain_model.py                  # Mesa Model: domain setup and scheduling
│   ├── element_agent_optimized.py       # ElementAgent: per-element cell/cytokine ODEs
│   ├── bone_healing_model_optimized.py  # ODE right-hand-side definitions
│   ├── endothelial_cell_agent.py        # Endothelial cell agent
│   ├── neighbor_cache_patch.py          # Neighbor-lookup optimization
│   └── utils.py                         # Mesh parsing utilities
├── output/
│   └── output.csv          # Simulation output (generated on run)
├── environment.yml          # Conda environment specification
└── README.md
```

---

## Requirements

- [Conda](https://docs.conda.io/en/latest/) (Miniconda or Anaconda)
- Python 3.13

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/ahedayat99/IMMUNO-BONE.git
cd IMMUNO-BONE
```

### 2. Create and activate the conda environment

```bash
conda env create -f environment.yml
conda activate ABM_env
```

This installs all required dependencies including:
- **Mesa 2.4** — agent-based modeling framework
- **NumPy / SciPy** — numerical computation
- **Pandas / Matplotlib** — data handling and visualization
- **NetworkX** — graph/neighbor operations
- **Pathos / Multiprocess** — parallel processing support

---

## Running the Model

All commands should be run from the **repository root directory**.

### Basic run (with calibrated parameters)

```bash
python scripts/simple_domain_modular_v1.py \
    --params_json data/input_params.json \
    --node data/node_elements.txt \
    --output_csv output/output.csv
```

### Verbose run (prints hourly state to console)

```bash
python scripts/simple_domain_modular_v1.py \
    --params_json data/input_params.json \
    --node data/node_elements.txt \
    --output_csv output/output.csv \
    --verbose
```


### Command-line arguments

| Argument | Description | Default |
|---|---|---|
| `--node` | Path to the FE mesh file (required) | — |
| `--params_json` | Path to JSON parameter file | Built-in defaults |
| `--output_csv` | Path for the output CSV file | `simulation_outputs.csv` |
| `--verbose` | Print hourly progress to stdout | Off |

---

## Input Files

### `data/node_elements.txt`

An Abaqus-format mesh file defining the 2D callus geometry. Contains node coordinates and element connectivity. The model reads this to construct the spatial domain over which agents are distributed.

### `data/input_params.json`

A JSON dictionary of calibrated kinetic and rate parameters governing cell recruitment, polarization, apoptosis, and cytokine production/degradation. If omitted, the model uses the built-in default parameters defined in `simple_domain_modular_v1.py`.


## Output

The simulation runs for **120 hours** and writes one row per hour to the output CSV file.

### Output columns

| Column | Description |
|---|---|
| `hour` | Simulation time (hours post-fracture) |
| `total_PMN` | Total neutrophil concentration |
| `total_M0` | Total resting macrophage concentration |
| `total_M1` | Total pro-inflammatory macrophage concentration |
| `total_M2` | Total anti-inflammatory macrophage concentration |
| `total_MSC` | Total mesenchymal stem cell concentration |
| `total_EC` | Total endothelial cell concentration |
| `total_c1` | Total pro-inflammatory cytokine (TNF-α / IL-1β proxy) |
| `total_c2` | Total anti-inflammatory cytokine (IL-10 proxy) |
| `total_c3` | Total MSC chemoattractant (PDGF proxy) |
| `total_c4` | Total angiogenic cytokine (VEGF proxy) |
| `total_debris` | Total cellular debris load |
| `num_agents` | Number of active element agents |

---

## Citation

If you use this model in your research, please cite:

> [Paper citation — to be added upon publication]

---
