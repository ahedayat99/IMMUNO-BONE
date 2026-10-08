#!/bin/bash -l
#$ -P bonesgrp
#$ -l h_rt=24:00:00
#$ -pe omp 32
#$ -l mem_per_core=2G
#$ -N immuno_sa
#$ -cwd
#$ -j y
#$ -o logs/
# OAT sensitivity analysis (53 parameters x 3 eps x 2 signs x 3 seeds) at the option-B base point.
source scc/common.sh
OPTION=${OPTION:-B}
python -u tools/oat_sa.py --base params/option_${OPTION}/sa_base.json --out out/option_${OPTION}/sa --workers $WORKERS
