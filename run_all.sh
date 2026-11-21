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

# ground truths of simulator/make_truth.py, one folder data/reference/S/ each: painted
# (mu and alpha painted) and, for each R of RATIOS, painted_rR (mu of painted, constant r = R)
RATIOS="1.2 2 2.8"
# standard deviations (ms) of the Gaussian noise on the LATs of painted: for each N a scenario
# painted_noiseN, ground truth of painted, LATs of Simulate with SIGMA = N
NOISES="0.05 0.1 0.3"
# the noise of the core sample in anisotropy: scenarios painted_rR_noiseCORE_NOISE
CORE_NOISE=0.05

TRUTHS="painted"
for r in $RATIOS; do TRUTHS="$TRUTHS painted_r$r"; done

# the scenarios of the grid (painted, then painted_noiseN): per-face estimator and Roney
GRID_SCENARIOS="painted"
for n in $NOISES; do GRID_SCENARIOS="$GRID_SCENARIOS painted_noise$n"; done
# all the scenarios of the per-face estimator: the grid, then the core sample (per-face only)
ALL_SCENARIOS="$GRID_SCENARIOS"
for r in $RATIOS; do ALL_SCENARIOS="$ALL_SCENARIOS painted_r${r}_noise$CORE_NOISE"; done

cmake -S estimator -B estimator/build -DCMAKE_BUILD_TYPE=Release
cmake --build estimator/build

# ground truths in format D; the six LAT maps (one per row of pacing_sites.csv) of painted without noise
for t in $TRUTHS; do
    python simulator/make_truth.py $t                             # -> data/reference/$t/truth.ply
done
./estimator/build/Simulate data/reference/painted/truth.ply simulator/stage1/pacing_sites.csv \
    data/reference/painted/lat.ply 0                              # -> data/reference/painted/lat.ply

# the noisy scenarios: a ground truth copied, the LATs with noise of standard deviation SIGMA ms
noisy() {                                                         # noisy TRUTH SCENARIO SIGMA
    mkdir -p data/reference/$2
    cp data/reference/$1/truth.ply data/reference/$2/truth.ply
    ./estimator/build/Simulate data/reference/$2/truth.ply simulator/stage1/pacing_sites.csv \
        data/reference/$2/lat.ply $3                              # -> data/reference/$2/lat.ply
}
for n in $NOISES; do noisy painted painted_noise$n $n; done
for r in $RATIOS; do noisy painted_r$r painted_r${r}_noise$CORE_NOISE $CORE_NOISE; done

rm -rf output                                     # no old results in the figures

# for each scenario all subsets of at least 3 maps with the per-face estimator; on the scenarios
# of the grid only, all triplets with method 2 of Roney et al. 2019; the estimators read only lat.ply
for s in $ALL_SCENARIOS; do
    python experiments/subsets.py $s              # -> output/$s/confs/*.conf
done
for s in $ALL_SCENARIOS; do
    ./estimator/build/Estimate output/$s/confs/*.conf   # -> output/$s/per_face/sub_*/estimate.{ply,csv}
done
python estimator/roney.py experiments/roney.conf $GRID_SCENARIOS  # -> output/SCENARIO/roney/sub_i-j-k/estimate.{ply,csv}

# errors against data/reference/SCENARIO/truth.ply, tables and figures
python analysis/analysis.py                       # -> output/figures/

# PLYs for MeshLab, all in output/mesh/: the cases are chosen in error_maps.py, which writes
# the errors of each case next to its estimate.ply, in error_maps.txt, creates the folder
# output/mesh/METHOD_SCENARIO_SUB/ and prints "MESH NUMBERS CASEDIR". Per case: the error
# PLY (log-Euclidean tensor error) and FiberView on the estimate in CASEDIR/, FiberView on the ground truth of its scenario
# in output/mesh/truth_SCENARIO/ (rewritten once per case of the scenario, same files).
python experiments/error_maps.py | while read -r mesh numbers case_dir; do
    ./estimator/build/ErrorMaps "$mesh" "$numbers" "$case_dir/error.ply"   # -> error.ply
    ./estimator/build/FiberView "$(dirname "$numbers")/estimate.ply" "$case_dir/estimate"  # -> estimate_{ratio,fibers}.ply
    scenario=$(basename "$(dirname "$mesh")")
    mkdir -p "output/mesh/truth_$scenario"
    ./estimator/build/FiberView "$(dirname "$mesh")/truth.ply" "output/mesh/truth_$scenario/truth"  # -> truth_{ratio,fibers}.ply
done
