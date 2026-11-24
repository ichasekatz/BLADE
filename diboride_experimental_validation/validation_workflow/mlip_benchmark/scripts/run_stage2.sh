#!/bin/bash
# Stage 2 (reduced): every model sequentially on the GPU, one process at a time. Resumable.
cd "$(dirname "$0")/.."
E=/home/ichasekatz/BLADE/mlip_envs
for m in mattersim mace-mpa sevennet mace-omat; do
  case $m in mace*) py=$E/mace/.venv/bin/python;; mattersim) py=$E/mattersim/.venv/bin/python;; sevennet) py=$E/sevennet/.venv/bin/python;; esac
  echo "$(date) start $m"
  $py scripts/m01_relax_runner.py $m --manifest stage2_manifest.csv 2>&1 | grep --line-buffered "^\["
  echo "$(date) done $m"
done
