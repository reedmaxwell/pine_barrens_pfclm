#!/bin/bash
# The spin-up sequence (East River docs/SPINUP_HOWTO.md, the Delaware notes), end to end:
#   1. one coupled year from the current wet state, to get the deck's own P - ET
#   2. the signed P - ET as a top-layer source
#   3. steady solve, seepage face on top (guard on, growth dt with Newton control)
#   4. steady solve, kinematic top: the streams form
#   5. coupled recursion, N cycles, each from the previous end state, scored each time
#
#   bash scripts/respin.sh <tag> <start pressure> <start CLM restart> [cycles]
set -e
TAG=$1; IC=$2; RST=$3; CYCLES=${4:-4}
ROOT=$(cd "$(dirname "$0")/.." && pwd); cd "$ROOT"
PF_COUPLED=/Users/reed/parflow/build_er_trial                         # PR 777 + source pftools
PF_STEADY=/Users/reed/Projects/TFG_smoother/parflow/build_etguard_adt   # guard + adaptive dt
PY=/Users/reed/miniforge3.1/envs/subsettools/bin/python
mkdir -p runs/${TAG}_states
echo "=== 1. coupled year ${TAG}_y1 from $IC"
PARFLOW_DIR=$PF_COUPLED $PF_COUPLED/py-env/bin/python scripts/run_pine_barrens.py --name ${TAG}_y1 --ic "$IC" --clm-restart "$RST" | grep -E 'successfully|failed'
PYTHONPATH=$PF_COUPLED/python $PY scripts/evaluate.py runs/${TAG}_y1 2>&1 | grep -vE 'Warning|c /=|hour '
$PY scripts/plot_et.py runs/${TAG}_y1 2>&1 | grep -E 'mask-mean|wrote'
cp runs/${TAG}_y1/${TAG}_y1.out.press.08784.pfb runs/${TAG}_states/y1_press.pfb; cp runs/${TAG}_y1/washita.rst.00000.0 runs/${TAG}_states/y1_clm.rst
echo "=== 2. P - ET"
$PY scripts/make_pme.py runs/${TAG}_y1 --out inputs/static/pme_${TAG}.pfb
echo "=== 3. steady, seepage face"
rm -rf runs/${TAG}_seep; PARFLOW_DIR=$PF_STEADY PYTHONPATH=$PF_STEADY/python $PY scripts/run_steady.py --ic runs/${TAG}_states/y1_press.pfb --name ${TAG}_seep --top seepage --pme inputs/static/pme_${TAG}.pfb | grep -E 'successfully|failed|top '
$PY scripts/check_steady.py runs/${TAG}_seep | tail -2
cp $(ls runs/${TAG}_seep/${TAG}_seep.out.press.*.pfb | sort | tail -1) runs/${TAG}_states/seep_end.pfb
echo "=== 4. steady, kinematic top"
rm -rf runs/${TAG}_overland; PARFLOW_DIR=$PF_STEADY PYTHONPATH=$PF_STEADY/python $PY scripts/run_steady.py --ic runs/${TAG}_states/seep_end.pfb --name ${TAG}_overland --top overland --pme inputs/static/pme_${TAG}.pfb | grep -E 'successfully|failed|top '
$PY scripts/check_steady.py runs/${TAG}_overland | tail -2
cp $(ls runs/${TAG}_overland/${TAG}_overland.out.press.*.pfb | sort | tail -1) runs/${TAG}_states/overland_end.pfb
echo "=== 5. coupled recursion, $CYCLES cycles"
$PY scripts/spinup.py --cycles $CYCLES --tag ${TAG}_rec --start runs/${TAG}_states/overland_end.pfb --clm-start runs/${TAG}_states/y1_clm.rst 2>&1 | grep -vE 'Warning|c /=|hour '
echo "=== respin $TAG done"
