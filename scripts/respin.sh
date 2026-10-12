#!/bin/zsh
# The spin-up sequence (East River docs/SPINUP_HOWTO.md, the Delaware notes), end to end:
#   1. one coupled year from the current wet state, to get the deck's own P - ET
#   2. the signed P - ET as a top-layer source
#   3. steady solve, seepage face on top (guard on, growth dt with Newton control)
#   4. steady solve, kinematic top: the streams form
#   5. coupled recursion, N cycles, each from the previous end state, scored each time
#
#   bash scripts/respin.sh <tag> <start pressure> <start CLM restart or none> [cycles]
#   RUN_ARGS="--geology-k-scale 100 --horizontal-only --flow-barrier $INPUTS/static/pf_flowbarrier_42m.pfb" for deck variants
#   (the steady solves take the same subsurface through STEADY_ARGS, see run_steady.py)
set -e
TAG=$1; IC=$2; RST=$3; CYCLES=${4:-4}; RUN_ARGS=${RUN_ARGS:-}; STEADY_ARGS=${STEADY_ARGS:-}
ROOT=$(cd "$(dirname "$0")/.." && pwd); cd "$ROOT"
DOM=${PB_DOMAIN:-v1}; export PB_DOMAIN=$DOM; RUNS=runs/$DOM; INPUTS=inputs/$DOM
PF_COUPLED=/Users/reed/parflow/build_er_trial                         # PR 777 + source pftools
PF_STEADY=/Users/reed/Projects/TFG_smoother/parflow/build_etguard_adt   # guard + adaptive dt
PY=/Users/reed/miniforge3.1/envs/subsettools/bin/python
mkdir -p $RUNS/${TAG}_states
START_STAGE=${START_STAGE:-1}                # START_STAGE=2 resumes after a finished year 1
if [ "$START_STAGE" -le 1 ]; then
echo "=== 1. coupled year ${TAG}_y1 from $IC"
if [ "$RST" = "none" ]; then RSTARG=""; else RSTARG="--clm-restart $RST"; fi     # "none": CLM cold start
PARFLOW_DIR=$PF_COUPLED $PF_COUPLED/py-env/bin/python scripts/run_pine_barrens.py --name ${TAG}_y1 --ic "$IC" ${=RSTARG} ${=RUN_ARGS} | grep -E 'successfully|failed'
cp $RUNS/${TAG}_y1/${TAG}_y1.out.press.08784.pfb $RUNS/${TAG}_states/y1_press.pfb; cp $RUNS/${TAG}_y1/washita.rst.00000.0 $RUNS/${TAG}_states/y1_clm.rst
PYTHONPATH=$PF_COUPLED/python $PY scripts/evaluate.py $RUNS/${TAG}_y1 2>&1 | grep -vE 'Warning|c /=|hour '
$PY scripts/plot_et.py $RUNS/${TAG}_y1 2>&1 | grep -E 'mask-mean|wrote' || true
fi
echo "=== 2. P - ET"
$PY scripts/make_pme.py $RUNS/${TAG}_y1 --out $INPUTS/static/pme_${TAG}.pfb
echo "=== 3. steady, seepage face"
rm -rf $RUNS/${TAG}_seep; PARFLOW_DIR=$PF_STEADY PYTHONPATH=$PF_STEADY/python $PY scripts/run_steady.py --ic $RUNS/${TAG}_states/y1_press.pfb --name ${TAG}_seep --top seepage --pme $INPUTS/static/pme_${TAG}.pfb ${=STEADY_ARGS} | grep -E 'successfully|failed|top '
$PY scripts/check_steady.py $RUNS/${TAG}_seep | tail -2
cp $(ls $RUNS/${TAG}_seep/${TAG}_seep.out.press.*.pfb | sort | tail -1) $RUNS/${TAG}_states/seep_end.pfb
echo "=== 4. steady, kinematic top"
rm -rf $RUNS/${TAG}_overland; PARFLOW_DIR=$PF_STEADY PYTHONPATH=$PF_STEADY/python $PY scripts/run_steady.py --ic $RUNS/${TAG}_states/seep_end.pfb --name ${TAG}_overland --top overland --pme $INPUTS/static/pme_${TAG}.pfb ${=STEADY_ARGS} | grep -E 'successfully|failed|top '
$PY scripts/check_steady.py $RUNS/${TAG}_overland | tail -2
cp $(ls $RUNS/${TAG}_overland/${TAG}_overland.out.press.*.pfb | sort | tail -1) $RUNS/${TAG}_states/overland_end.pfb
echo "=== 5. coupled recursion, $CYCLES cycles"
$PY scripts/spinup.py --cycles $CYCLES --tag ${TAG}_rec --start $RUNS/${TAG}_states/overland_end.pfb --clm-start $RUNS/${TAG}_states/y1_clm.rst --run-args "$RUN_ARGS" 2>&1 | grep -vE 'Warning|c /=|hour '
echo "=== respin $TAG done"
