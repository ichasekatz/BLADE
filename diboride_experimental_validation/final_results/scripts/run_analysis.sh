#!/bin/bash
# Full analysis chain on the cached MLIP energies (no relaxations):
# z01 ordered hull / CE -> x01 six-system thermodynamics -> x02 sensitivity -> d01 deliverables -> v01 validation
# -> g01 figures. Logs to final_results/logs/analysis_<step>.log.
set -e
B=/home/ichasekatz/BLADE/BLADE
S=$B/diboride_experimental_validation/final_results/scripts
L=$B/diboride_experimental_validation/final_results/logs
export MPLCONFIGDIR=$HOME/.cache/matplotlib
cd $B
for step in z01_vacancy_ordering_ce_mc x01_six_system_thermo x02_six_system_sensitivity d01_deliverables v01_validation g01_figures; do
  echo "$(date) start $step" >> $L/analysis_chain.log
  .venv/bin/python $S/$step.py > $L/analysis_$step.log 2>&1
  echo "$(date) done $step" >> $L/analysis_chain.log
done
echo "$(date) ALL DONE" >> $L/analysis_chain.log
