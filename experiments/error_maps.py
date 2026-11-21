#!/usr/bin/env python3
"""
Error maps in PLY for MeshLab, from the runs of the grid. The error, the
log-Euclidean distance between estimated and true tensor, comes from load() of
analysis/analysis.py, against the ground truth.
For each case writes next to its estimate.ply the file error_maps.txt, only
numbers: one row per face, in face_id order, "estimable error", 0 0 on the
non-estimable faces. Then creates the folder of the case,
CASEDIR = output/mesh/METHOD_SCENARIO_SUB, and prints a line
"MESH NUMBERS CASEDIR": MESH is data/reference/SCENARIO/lat.ply, NUMBERS the
error_maps.txt. run_all.sh passes them to ErrorMaps with OUT CASEDIR/error.ply,
which writes the error as the quality of the faces. The cases, for each
scenario of GRID (painted at every LAT noise level, from analysis.py):
    per-face, all maps
    per-face and Roney on one triplet, the same at every level so that in
    MeshLab the noise levels compare on the same maps: the triplet with the
    largest paired area (faces estimable by both methods, area-weighted) on
    painted without noise, the first scenario of GRID

Usage, from the repository root:
    python experiments/error_maps.py | while read -r mesh numbers case_dir; do
        ./estimator/build/ErrorMaps "$mesh" "$numbers" "$case_dir/error.ply"
    done
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]                # repository root
sys.path.insert(0, str(ROOT / "analysis"))
from analysis import GRID, load, subs                     # the errors of a run, the scenarios


def line(method, scenario, sub, df):
    """Writes error_maps.txt of the case, creates its folder in output/mesh/ and prints its line. Numbers with repr, which reads back
    exactly as a double; no "nan", which >> in C++ does not read."""
    rows = []
    for e, err in zip(df["estimable"].tolist(), df["error"].tolist()):
        rows.append(f"1 {err!r}" if e else "0 0")
    case = f"output/{scenario}/{method}/{sub}"
    (ROOT / case / "error_maps.txt").write_text("\n".join(rows) + "\n")
    case_dir = f"output/mesh/{method}_{scenario}_{sub}"
    (ROOT / case_dir).mkdir(parents=True, exist_ok=True)
    print(f"data/reference/{scenario}/lat.ply {case}/error_maps.txt {case_dir}")


def paired_area(scenario, sub):
    """Area (mm^2) of the faces estimable by both methods with the triplet sub."""
    per_face, roney = load(scenario, "per_face", sub), load(scenario, "roney", sub)
    return roney["area"][roney["estimable"] & per_face["estimable"]].sum()


def main():
    clean = GRID[0][0]                                    # painted without noise: chooses the triplet
    triplet = max(subs(clean, "roney"), key=lambda sub: paired_area(clean, sub))
    for scenario, _ in GRID:
        sub = subs(scenario, "per_face")[-1]              # all maps
        line("per_face", scenario, sub, load(scenario, "per_face", sub))
        line("per_face", scenario, triplet, load(scenario, "per_face", triplet))
        line("roney", scenario, triplet, load(scenario, "roney", triplet))


if __name__ == "__main__":
    main()
