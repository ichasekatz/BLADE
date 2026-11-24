#!/bin/bash
# MACE-MPA, SevenNet and MACE-OMAT Stage 4 in three tiers so the six-system thermodynamics can run early:
#   1. essential  (58): hosts, x=0.5 vacancy-environment cells, end-member vacancies, spread multi-vacancy cells,
#                       clay1 ordered cells, competitor references
#   2. priority  (102): + x=0.25/0.75 vacancy cells, coupled-layer cells, clusters/planes, other configurations
#   3. full      (548): + enumerated ordered TaB2-x / NbB2-y cells
# Resumable: the runner skips structures with an existing result.json.
B=/home/ichasekatz/BLADE/BLADE/diboride_experimental_validation
SD=$B/final_results/supplementary_data
LG=$B/final_results/logs
E=/home/ichasekatz/BLADE/mlip_envs
export MPLCONFIGDIR=$HOME/.cache/matplotlib
cd $B/validation_workflow/mlip_benchmark
tiers() {
  local m=$1 py=$2
  for t in essential priority; do
    $py scripts/m01_relax_runner.py $m --manifest $SD/stage4_manifest_$t.csv
    echo "$(date) tier $t done $m"
  done
  $py scripts/m01_relax_runner.py $m --manifest $SD/stage4_manifest.csv
  echo "$(date) all done $m"
}
tiers mace-mpa $E/mace/.venv/bin/python >> $LG/stage4_macempa.log 2>&1 &
tiers sevennet $E/sevennet/.venv/bin/python >> $LG/stage4_sevennet.log 2>&1 &
tiers mace-omat $E/mace/.venv/bin/python >> $LG/stage4_maceomat_early.log 2>&1 &
wait
