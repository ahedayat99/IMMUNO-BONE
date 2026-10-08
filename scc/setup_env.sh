#!/bin/bash -l
# One-time environment setup on SCC (run interactively on a login node):  bash scc/setup_env.sh
module load miniconda
conda create -y -n immuno_v2 python=3.12 numpy scipy pandas
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate immuno_v2
pip install "mesa==2.4.0"
python -c "import mesa, numpy, scipy; print('mesa', mesa.__version__, 'numpy', numpy.__version__, 'scipy', scipy.__version__)"
