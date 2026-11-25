#!/usr/bin/env python3
"""
The cases of the PLYs for MeshLab, all in output/mesh/, in three folders. The error,
the log-Euclidean distance between estimated and true tensor, comes from load() of
analysis/analysis.py, against the ground truth.

    roney_comparison/  for each scenario of GRID (painted at every LAT noise level, from
                       analysis.py), per-face and Roney on one triplet, the same at every
                       level so that the levels compare on the same maps: the triplet with
                       the largest paired area (faces estimable by both methods,
                       area-weighted) on painted without noise, the first scenario of GRID
                       -> METHOD_SCENARIO_SUB_error.ply
    worst_subsets/     on the scenario of GRID with sigma = WORST_NOISE, per-face: for each
                       k = 3..6 the subset with k maps with the highest mean error,
                       area-weighted on its estimable faces
                       -> per_face_SCENARIO_kK_SUB_error.ply
    ellipses/          the fields as ellipses: the ground truth of painted, and on the
                       scenario of worst_subsets the per-face estimate of the worst triplet
                       and of all the maps
                       -> truth_painted_ellipses.ply, per_face_SCENARIO_kK_SUB_ellipses.ply

For each error map writes next to its estimate.ply the file error_maps.txt, only numbers:
one row per face, in face_id order, "estimable error", 0 0 on the non-estimable faces.
Creates the three folders and prints one line per case, the kind first:
    error MESH NUMBERS OUT      MESH = data/reference/SCENARIO/lat.ply, NUMBERS its error_maps.txt
    ellipses D OUT              D = a truth.ply or estimate.ply (format D)
with OUT the PLY to write. run_all.sh passes them to ErrorMaps and FiberView.

Usage, from the repository root:
    python experiments/error_maps.py | while read -r kind a b c; do
        if [ "$kind" = error ]; then ./estimator/build/ErrorMaps "$a" "$b" "$c"
        elif [ "$kind" = ellipses ]; then ./estimator/build/FiberView "$a" "$b"; fi
    done
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]                # repository root
sys.path.insert(0, str(ROOT / "analysis"))
from analysis import GRID, all_subsets, load, subs        # the errors of the runs, the scenarios

MESH = "output/mesh"
WORST_NOISE = 0.05                                        # sigma (ms) of the scenario of worst_subsets
KS = range(3, 7)                                          # maps of the subsets of worst_subsets


def error_line(folder, method, scenario, sub, name):
    """Writes error_maps.txt of the case and prints its line, OUT = MESH/FOLDER/NAME_error.ply. Numbers
    with repr, which reads back exactly as a double; no "nan", which >> in C++ does not read."""
    df = load(scenario, method, sub)
    rows = []
    for e, err in zip(df["estimable"].tolist(), df["error"].tolist()):
        rows.append(f"1 {err!r}" if e else "0 0")
    case = f"output/{scenario}/{method}/{sub}"
    (ROOT / case / "error_maps.txt").write_text("\n".join(rows) + "\n")
    print(f"error data/reference/{scenario}/lat.ply {case}/error_maps.txt {MESH}/{folder}/{name}_error.ply")


def paired_area(scenario, sub):
    """Area (mm^2) of the faces estimable by both methods with the triplet sub."""
    per_face, roney = load(scenario, "per_face", sub), load(scenario, "roney", sub)
    return roney["area"][roney["estimable"] & per_face["estimable"]].sum()


def mean_errors(scenario, k):
    """Per-face, the subsets with k maps of the scenario: mean error of each, sum(area * error) /
    sum(area) on its estimable faces, as a Series indexed by sub."""
    d = all_subsets(scenario)
    est = d[(d["k"] == k) & d["estimable"]]
    return (est["area"] * est["error"]).groupby(est["sub"]).sum() / est["area"].groupby(est["sub"]).sum()


def main():
    for folder in ("roney_comparison", "worst_subsets", "ellipses"):
        (ROOT / MESH / folder).mkdir(parents=True, exist_ok=True)

    clean = GRID[0][0]                                    # painted without noise: chooses the triplet
    triplet = max(subs(clean, "roney"), key=lambda sub: paired_area(clean, sub))
    for scenario, _ in GRID:
        for method in ("per_face", "roney"):
            error_line("roney_comparison", method, scenario, triplet, f"{method}_{scenario}_{triplet}")

    noisy = next(sc for sc, n in GRID if n == WORST_NOISE)
    worst = {k: mean_errors(noisy, k).idxmax() for k in KS}
    for k, sub in worst.items():
        error_line("worst_subsets", "per_face", noisy, sub, f"per_face_{noisy}_k{k}_{sub}")

    print(f"ellipses data/reference/{clean}/truth.ply {MESH}/ellipses/truth_{clean}_ellipses.ply")
    for k in (KS[0], KS[-1]):                             # the worst triplet, all the maps
        print(f"ellipses output/{noisy}/per_face/{worst[k]}/estimate.ply "
              f"{MESH}/ellipses/per_face_{noisy}_k{k}_{worst[k]}_ellipses.ply")


if __name__ == "__main__":
    main()
