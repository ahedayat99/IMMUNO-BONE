# Sourced by the job scripts. Edit ROOT if the repository is uploaded elsewhere.
ROOT="${SGE_O_WORKDIR:-$PWD}"   # the directory the job was submitted from (repository root)
module load miniconda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate immuno_v2
export PYTHONUNBUFFERED=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$ROOT" || { echo "cannot cd to $ROOT"; exit 1; }
python -c "import mesa, numpy, scipy; assert mesa.__version__.startswith('2.4'), 'need mesa 2.4.x, found ' + mesa.__version__; print('env ok: mesa', mesa.__version__)" || exit 1
WORKERS=$(( ${NSLOTS:-8} - 1 ))
