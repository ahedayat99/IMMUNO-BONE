# IMMUNO-BONE: Bone Healing Agent-Based Model (ABM)

This repository contains the full implementation of the IMMUNO-BONE model — a hybrid agent-based / continuum simulation of the early inflammatory and regenerative phases of murine bone fracture healing. The model captures spatial interactions among key cell populations (neutrophils, macrophages M0/M1/M2, mesenchymal stromal/progenitor cells, and endothelial cells) and pro- and anti-inflammatory cytokines (TNF-α, IL-10, TGF-β, VEGF), distributed across a finite-element mesh of the osteotomy gap and periosteal callus.

**Version 2.0.0** is the revised model accompanying the revised manuscript (npj Systems Biology and Applications). Version 1.0.x (the originally submitted model) remains available under the `v1.0.0` / `V1.0.1` tags and on Zenodo.

### What changed in v2.0.0

- **Molecular fields:** cytokines are transported between elements by a mass-conserving finite-volume diffusion operator (`scripts/mesh.py`); diffusion coefficients in mm²/h.
- **Reaction step:** cells, cytokines and debris in each element are integrated together as one coupled ODE system over each 1-h step (`scripts/reaction.py`, `scripts/parameters.py`).
- **Agent layer:** phenotype (composition)-weighted chemotactic migration of cell parcels, per-element capacity, and removed parcels are deleted cleanly from element occupancy.
- **Equations:** the fibrogenic term F1 enters with a negative sign, and VEGF has a sink proportional to its concentration (see the manuscript Methods).
- **Initial conditions:** progenitor (MSC) layout switch `init_MSC_layout` (0 = gap, 1 = periosteal, 2 = half periosteal / half gap; the model of record uses 2).
- **Recalibration:** joint calibration to the immunofluorescence macrophage ratios and the dPCR cytokine fold-changes, with literature-derived parameters kept inside their published ranges → `params/model_of_record.json`.
- **Reproducibility:** explicit random seeds (`--seed`).

---

## Repository Structure

```
.
├── data/
│   └── node_elements.txt          # FE mesh geometry (nodes + elements) of the callus
├── params/
│   └── model_of_record.json       # Calibrated parameter set used for all reported results
├── scripts/
│   ├── simple_domain_modular_v1.py      # Main entry point — run this
│   ├── domain_model.py                  # Mesa Model: domain setup, reaction/transport/migration scheduling
│   ├── element_agent_optimized.py       # ElementAgent: cell parcels and budding rules
│   ├── reaction.py                      # Element-local coupled ODE (cells + cytokines + debris)
│   ├── bone_healing_model_optimized.py  # ODE right-hand-side terms
│   ├── mesh.py                          # Mesh adjacency, element areas, finite-volume diffusion
│   ├── zones.py                         # Anatomical zones matched to the histological ROIs
│   ├── parameters.py                    # Time-step and unit conventions
│   ├── endothelial_cell_agent.py        # Endothelial cell agent
│   ├── neighbor_cache_patch.py          # Neighbor-lookup optimization
│   └── utils.py                         # Mesh parsing utilities
├── environment.yml                # Conda environment specification
├── LICENSE                        # MIT License
├── CITATION.cff
└── README.md
```

---

## Requirements

- [Conda](https://docs.conda.io/en/latest/) (Miniconda or Anaconda)
- Python 3.13, Mesa 2.4

---

## Installation

```bash
git clone https://github.com/ahedayat99/IMMUNO-BONE.git
cd IMMUNO-BONE
conda env create -f environment.yml
conda activate ABM_env
```

This installs all required dependencies including Mesa 2.4, NumPy / SciPy, Pandas / Matplotlib, NetworkX and Pathos / Multiprocess.

---

## Running the Model

All commands should be run from the **repository root directory**. One 120-h run takes roughly 10–40 s on a laptop.

### Model of record (as reported in the revised manuscript)

```bash
python scripts/simple_domain_modular_v1.py \
    --params_json params/model_of_record.json \
    --node data/node_elements.txt \
    --output_csv output/output.csv \
    --seed 11000
```

Add `--verbose` to print the hourly state. The reported results are means over the 10 seeds 11000, 12000, …, 20000.

### Command-line arguments

| Argument | Description | Default |
|---|---|---|
| `--node` | Path to the FE mesh file (required) | — |
| `--params_json` | Path to JSON parameter file | Built-in defaults |
| `--output_csv` | Path for the output CSV file | `simulation_outputs.csv` |
| `--seed` | Random seed; the same seed reproduces a run exactly | Unseeded |
| `--verbose` | Print hourly progress to stdout | Off |

---

## Input Files

### `data/node_elements.txt`

An Abaqus-format mesh file defining the 2D callus geometry (3,242 eight-node quadrilaterals): the osteotomy gap and the periosteal callus in front of the far cortex. Coordinates are in mm, origin at the centre of the gap; x is perpendicular to the fixation plate and y runs along the bone axis.

### `params/model_of_record.json`

The calibrated parameter set: kinetic and rate parameters governing cell recruitment, polarization, apoptosis, and cytokine production/degradation/diffusion (rates per hour). Parameters missing from a file take the built-in defaults defined in `simple_domain_modular_v1.py`.

---

## Output

The simulation runs for **120 hours** and writes one row per hour to the output CSV file.

| Column | Description |
|---|---|
| `hour` | Simulation time (hours post-fracture) |
| `total_PMN` | Total neutrophils |
| `total_M0` | Total resting macrophages |
| `total_M1` | Total pro-inflammatory macrophages |
| `total_M2` | Total anti-inflammatory macrophages |
| `total_MSC` | Total mesenchymal stromal/progenitor cells |
| `total_EC` | Total endothelial cells |
| `total_c1` | Total pro-inflammatory cytokine (TNF-α) |
| `total_c2` | Total anti-inflammatory cytokine (IL-10) |
| `total_c3` | Total TGF-β |
| `total_c4` | Total angiogenic cytokine (VEGF) |
| `total_debris` | Total cellular debris load |
| `num_agents` | Number of active element agents (cell parcels) |

---

## Citation

If you use this model in your research, please cite the software (all versions, Zenodo concept DOI):
DOI: 10.5281/zenodo.21567254

> Article citation to be added upon publication.

---

## License

Released under the [MIT License](LICENSE). Copyright (c) 2026 Ahmad Hedayatzadeh Razavi and the Musculoskeletal Translational Innovation Initiative (MTII).
