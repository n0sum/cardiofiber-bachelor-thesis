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
# the number of maps of the core sample: only the subsets with CORE_K maps; equal to CORE_K in analysis.py
CORE_K=4

TRUTHS="painted"
for r in $RATIOS; do TRUTHS="$TRUTHS painted_r$r"; done

# the scenarios of the grid (painted, then painted_noiseN): per-face estimator and Roney
GRID_SCENARIOS="painted"
for n in $NOISES; do GRID_SCENARIOS="$GRID_SCENARIOS painted_noise$n"; done
# the scenarios of the core sample in anisotropy (painted_rR_noiseCORE_NOISE): per-face only, k = CORE_K
CORE_SCENARIOS=""
for r in $RATIOS; do CORE_SCENARIOS="$CORE_SCENARIOS painted_r${r}_noise$CORE_NOISE"; done
# all the scenarios of the per-face estimator: the grid, then the core sample
ALL_SCENARIOS="$GRID_SCENARIOS$CORE_SCENARIOS"

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

# per-face estimator: on the scenarios of the grid all subsets of at least 3 maps (42), on those
# of the core sample only the subsets with CORE_K maps (15); on the scenarios of the grid only,
# all triplets with method 2 of Roney et al. 2019; the estimators read only lat.ply
for s in $GRID_SCENARIOS; do
    python experiments/subsets.py $s              # -> output/$s/confs/*.conf, 42
done
for s in $CORE_SCENARIOS; do
    python experiments/subsets.py $s $CORE_K      # -> output/$s/confs/*.conf, 15
done
for s in $ALL_SCENARIOS; do
    ./estimator/build/Estimate output/$s/confs/*.conf   # -> output/$s/per_face/sub_*/estimate.{ply,csv}
done
python estimator/roney.py experiments/roney.conf $GRID_SCENARIOS  # -> output/SCENARIO/roney/sub_i-j-k/estimate.{ply,csv}

# errors against data/reference/SCENARIO/truth.ply, tables and figures
python analysis/analysis.py                       # -> output/figures/

# PLYs for MeshLab, all in output/mesh/{roney_comparison,worst_subsets,ellipses}/: the cases are
# chosen in error_maps.py, which writes the errors of each error map next to its estimate.ply, in
# error_maps.txt, creates the three folders and prints one line per case, the kind first:
# "error MESH NUMBERS OUT" -> ErrorMaps (log-Euclidean error, one RGB color per face), "ellipses D OUT"
# -> FiberView (the field of D as ellipses on the gray mesh).
python experiments/error_maps.py | while read -r kind a b c; do
    if [ "$kind" = error ]; then
        ./estimator/build/ErrorMaps "$a" "$b" "$c"                # -> OUT, error map
    elif [ "$kind" = ellipses ]; then
        ./estimator/build/FiberView "$a" "$b"                     # -> OUT, ellipses
    fi
done
