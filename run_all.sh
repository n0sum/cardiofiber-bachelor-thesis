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

python simulator/stage1/opencarp/mm_to_um.py      # -> data/reference/LA_um.*

# scenarios aniso_zone_N: N random zones (seed 0), one of the 25 CV combinations of CV_CSV per zone;
# d = sqrt(4 A / (pi N)) ~ 123, 44, 15, 5.0 mm, equal steps on a log scale; N = 1: uniform CVs,
# N = 600: zones as large as Roney's disc (5 mm), which assumes uniform CVs within it
ZONES="1 8 70 600"
CV_CSV=experiments/region_cv/aniso_random.csv
for n in $ZONES; do
    s=aniso_zone_$n; mesh=data/reference/$s/LA_um
    python simulator/stage1/opencarp/random_tags.py $s $CV_CSV 0 $n   # -> $mesh.*
    python simulator/stage1/opencarp/run_sims.py $s $CV_CSV $mesh     # -> data/reference/sims/$s/run_*/act.igb
    python simulator/convert/toPly.py $s $CV_CSV $mesh                # -> data/reference/LA_{LAT,GT}_$s.ply
done

cmake -S estimator -B estimator/build -DCMAKE_BUILD_TYPE=Release
cmake --build estimator/build

rm -rf output                                     # no old results in the figures

# grid: for each scenario all subsets of at least 3 maps with the per-face estimator,
# and all triplets with method 2 of Roney et al. 2019
for n in $ZONES; do
    python experiments/subsets.py aniso_zone_$n   # -> output/aniso_zone_$n/confs/*.conf
done
./estimator/build/Estimate output/aniso_zone_*/confs/*.conf   # -> output/SCENARIO/per_face/sub_*/errors.csv
# -> output/SCENARIO/roney/sub_i-j-k/errors.csv
python estimator/roney.py experiments/roney.conf $(for n in $ZONES; do echo aniso_zone_$n; done)
python analysis/analysis.py                       # -> output/figures/

# error maps in PLY: per-face with all maps on aniso_zone_1 and aniso_zone_600, and the two
# methods on the triplet of aniso_zone_600 with the largest paired area; three PLYs per case, one per error
# error_maps.py writes the errors of each case next to its errors.csv, in error_maps.txt
python experiments/error_maps.py | while read -r mesh numbers prefix; do
    ./estimator/build/ErrorMaps "$mesh" "$numbers" "$prefix"   # -> output/mesh/*_{angle,cv_l,ratio}.ply
done
