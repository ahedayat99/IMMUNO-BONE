#!/bin/bash
# Sensitivity analysis, then the three calibration ranges in parallel once it finishes.
#   bash scc/submit_pipeline.sh            (option B)
#   OPTION=A bash scc/submit_pipeline.sh   (option A; needs params/option_A/sa_base.json)
cd "$(dirname "$0")/.." && mkdir -p logs
OPTION=${OPTION:-B}
SA=$(qsub -terse -v OPTION=$OPTION scc/submit_sa.sh)
echo "sensitivity job $SA"
for R in small medium large; do
  J=$(qsub -terse -hold_jid $SA -v OPTION=$OPTION,RANGE=$R scc/submit_calib.sh)
  echo "calibration $R job $J (waits for $SA)"
done
