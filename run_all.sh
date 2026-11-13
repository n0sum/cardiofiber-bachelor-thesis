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
# standard deviations (ms) of the Gaussian noise on the LATs: for each scenario S and each N
# a scenario S_noiseN, same ground truth, LATs of Simulate with SIGMA = N
NOISES="0.1 0.2"

# all the scenarios of the grid: first the clean ones, then those of each noise level
ALL_SCENARIOS="$SCENARIOS"
for n in $NOISES; do
    for s in $SCENARIOS; do
        ALL_SCENARIOS="$ALL_SCENARIOS ${s}_noise$n"
    done
done

cmake -S estimator -B estimator/build -DCMAKE_BUILD_TYPE=Release
cmake --build estimator/build

# ground truth in format D and the six LAT maps, one per row of pacing_sites.csv, without noise
for s in $SCENARIOS; do
    python simulator/make_truth.py $s                             # -> data/reference/$s/truth.ply
    ./estimator/build/Simulate data/reference/$s/truth.ply simulator/stage1/pacing_sites.csv \
        data/reference/$s/lat.ply 0                               # -> data/reference/$s/lat.ply
done
# the noisy scenarios: the ground truth copied, the LATs with noise of standard deviation N ms
for n in $NOISES; do
    for s in $SCENARIOS; do
        mkdir -p data/reference/${s}_noise$n
        cp data/reference/$s/truth.ply data/reference/${s}_noise$n/truth.ply
        ./estimator/build/Simulate data/reference/${s}_noise$n/truth.ply simulator/stage1/pacing_sites.csv \
            data/reference/${s}_noise$n/lat.ply $n                # -> data/reference/${s}_noise$n/lat.ply
    done
done

rm -rf output                                     # no old results in the figures

# grid: for each scenario all subsets of at least 3 maps with the per-face estimator,
# and all triplets with method 2 of Roney et al. 2019; the estimators read only lat.ply
for s in $ALL_SCENARIOS; do
    python experiments/subsets.py $s              # -> output/$s/confs/*.conf
done
for s in $ALL_SCENARIOS; do
    ./estimator/build/Estimate output/$s/confs/*.conf   # -> output/$s/per_face/sub_*/estimate.{ply,csv}
done
python estimator/roney.py experiments/roney.conf $ALL_SCENARIOS   # -> output/SCENARIO/roney/sub_i-j-k/estimate.{ply,csv}

# errors against data/reference/SCENARIO/truth.ply, statistics and figures: the clean run, then
# one run per noise level on the scenarios S_noiseN
python analysis/analysis.py                       # -> output/figures/
for n in $NOISES; do
    python analysis/analysis.py $n                # -> output/figures_noise$n/
done

# PLYs for MeshLab, all in output/mesh/: the cases are chosen in error_maps.py, which writes
# the errors of each case next to its estimate.ply, in error_maps.txt, creates the folder
# output/mesh/METHOD_SCENARIO_SUB/ and prints "MESH NUMBERS CASEDIR". Per case: the three error
# PLYs and FiberView on the estimate in CASEDIR/, FiberView on the ground truth of its scenario
# in output/mesh/truth_SCENARIO/ (for painted rewritten once per case, same files).
# Once without noise and once per noise level (error_maps.py N, scenarios S_noiseN).
for n in "" $NOISES; do
    python experiments/error_maps.py $n | while read -r mesh numbers case_dir; do
        ./estimator/build/ErrorMaps "$mesh" "$numbers" "$case_dir/error"   # -> error_{angle,mean_cv,ratio}.ply
        ./estimator/build/FiberView "$(dirname "$numbers")/estimate.ply" "$case_dir/estimate"  # -> estimate_{ratio,fibers}.ply
        scenario=$(basename "$(dirname "$mesh")")
        mkdir -p "output/mesh/truth_$scenario"
        ./estimator/build/FiberView "$(dirname "$mesh")/truth.ply" "output/mesh/truth_$scenario/truth"  # -> truth_{ratio,fibers}.ply
    done
done
