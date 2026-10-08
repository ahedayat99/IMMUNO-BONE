#!/bin/bash -l
#$ -P bonesgrp
#$ -l h_rt=48:00:00
#$ -pe omp 32
#$ -l mem_per_core=2G
#$ -N immuno_calib
#$ -cwd
#$ -j y
#$ -o logs/
# LHS + differential evolution for one search range. Submit with: qsub -v RANGE=medium,OPTION=B scc/submit_calib.sh
source scc/common.sh
OPTION=${OPTION:-B}; RANGE=${RANGE:-medium}
python -u tools/calibrate.py --base params/option_${OPTION}/sa_base.json \
    --subset out/option_${OPTION}/sa/calibration_subset.json \
    --range $RANGE --out out/option_${OPTION}/calib_${RANGE} --workers $WORKERS
