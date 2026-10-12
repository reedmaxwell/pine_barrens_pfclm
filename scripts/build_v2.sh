#!/bin/zsh
# Build domain v2 from HydroData, step by step (docs/DOMAIN_REDESIGN.md), then pull the observations.
set -e
cd "$(dirname "$0")/.."
export PB_DOMAIN=v2
PY=/Users/reed/miniforge3.1/envs/subsettools/bin/python
$PY scripts/build_domain.py --step delineate --force
for st in statics slopes clm urban mannings barrier_geology sites forcing; do
  $PY scripts/build_domain.py --step $st
done
$PY scripts/get_obs.py
echo "=== v2 build done"
