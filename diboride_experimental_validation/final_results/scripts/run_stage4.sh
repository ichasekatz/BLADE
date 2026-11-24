#!/bin/bash
# Stage 4 MLIP relaxations (all five models; two GPU processes at a time). Resumable; skips finished structures.
MB=/home/ichasekatz/BLADE/BLADE/diboride_experimental_validation/validation_workflow/mlip_benchmark
MAN=/home/ichasekatz/BLADE/BLADE/diboride_experimental_validation/final_results/supplementary_data/stage4_manifest.csv
E=/home/ichasekatz/BLADE/mlip_envs
cd $MB
run() { local m=$1 py=$2; echo "$(date) start $m"; $py scripts/m01_relax_runner.py $m --manifest $MAN 2>&1 | grep --line-buffered "^\[" ; echo "$(date) done $m"; }
export MPLCONFIGDIR=$HOME/.cache/matplotlib
run orb /home/ichasekatz/BLADE/BLADE/.venv/bin/python & run mattersim $E/mattersim/.venv/bin/python & wait
run mace-mpa $E/mace/.venv/bin/python & run sevennet $E/sevennet/.venv/bin/python & wait
run mace-omat $E/mace/.venv/bin/python
echo "$(date) ALL DONE"
