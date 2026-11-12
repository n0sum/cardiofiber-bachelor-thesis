#!/usr/bin/env bash
#
# Rebuilds data/reference/ and output/ from the AugmentA output
# (data/augmenta/, in the repository). Must be run with the conda environment cardiofiber active.
#
# Regenerating the AugmentA output from the SSM is optional and stays manual:
# see README.
#
set -euo pipefail
cd "$(dirname "$0")"

# scenarios: ground truth of simulator/make_truth.py, one folder data/reference/SCENARIO/ each
SCENARIOS="uniform painted r1.6 r2.4"

cmake -S estimator -B estimator/build -DCMAKE_BUILD_TYPE=Release
cmake --build estimator/build

# ground truth in format D and the six LAT maps, one per row of pacing_sites.csv
for s in $SCENARIOS; do
    python simulator/make_truth.py $s                             # -> data/reference/$s/truth.ply
    ./estimator/build/Simulate data/reference/$s/truth.ply simulator/stage1/pacing_sites.csv \
        data/reference/$s/lat.ply                                 # -> data/reference/$s/lat.ply
done

rm -rf output                                     # no old results in the figures

# grid: for each scenario all subsets of at least 3 maps with the per-face estimator,
# and all triplets with method 2 of Roney et al. 2019; the estimators read only lat.ply
for s in $SCENARIOS; do
    python experiments/subsets.py $s              # -> output/$s/confs/*.conf
done
for s in $SCENARIOS; do
    ./estimator/build/Estimate output/$s/confs/*.conf   # -> output/$s/per_face/sub_*/estimate.{ply,csv}
done
python estimator/roney.py experiments/roney.conf $SCENARIOS   # -> output/SCENARIO/roney/sub_i-j-k/estimate.{ply,csv}

# errors against data/reference/SCENARIO/truth.ply, statistics and figures
python analysis/analysis.py                       # -> output/figures/

# error maps in PLY: the cases are chosen in error_maps.py; three PLYs per case, one per error.
# error_maps.py writes the errors of each case next to its estimate.ply, in error_maps.txt
# For each case FiberView also writes the estimate and the ground truth of its scenario as PLYs
# for MeshLab, next to their D file; for painted the ground truth is rewritten once per case.
python experiments/error_maps.py | while read -r mesh numbers prefix; do
    ./estimator/build/ErrorMaps "$mesh" "$numbers" "$prefix"   # -> output/mesh/*_{angle,cv_l,ratio}.ply
    case_dir=$(dirname "$numbers")
    truth_dir=$(dirname "$mesh")
    ./estimator/build/FiberView "$case_dir/estimate.ply" "$case_dir/estimate"  # -> estimate_{ratio,fibers}.ply
    ./estimator/build/FiberView "$truth_dir/truth.ply" "$truth_dir/truth"      # -> truth_{ratio,fibers}.ply
done
