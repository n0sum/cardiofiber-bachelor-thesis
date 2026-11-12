#!/usr/bin/env python3
"""
Error maps in PLY for MeshLab, from the runs of the grid. The errors come from
load() of analysis/analysis.py, against the ground truth.
For each case writes next to its estimate.ply the file error_maps.txt, only
numbers: one row per face, in face_id order, "estimable angle cv_l ratio",
0 0 0 0 on the non-estimable faces. Then prints a line "MESH NUMBERS PREFIX",
to be passed to ErrorMaps, which writes PREFIX_{angle,cv_l,ratio}.ply.
MESH is data/reference/SCENARIO/lat.ply, NUMBERS the error_maps.txt, PREFIX
output/mesh/METHOD_SCENARIO_SUB. The cases:
    per-face, all maps, on each scenario of SCENARIOS_ALL_MAPS
    per-face and Roney on the triplet of SCENARIO_TRIPLET with the largest
    paired area: faces estimable by both methods, area-weighted

Usage, from the repository root:
    python experiments/error_maps.py | while read -r mesh numbers prefix; do
        ./estimator/build/ErrorMaps "$mesh" "$numbers" "$prefix"
    done
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]                # repository root
sys.path.insert(0, str(ROOT / "analysis"))
from analysis import load, subs                           # the errors of a run

SCENARIOS_ALL_MAPS = ["uniform", "painted"]
SCENARIO_TRIPLET = "painted"


def line(method, scenario, sub, df):
    """Writes error_maps.txt of the case and prints its line. Numbers with repr, which reads back
    exactly as a double; no "nan", which >> in C++ does not read."""
    rows = []
    for e, a, l, r in zip(df["estimable"].tolist(), df["angle"].tolist(), df["cv_l"].tolist(),
                          df["ratio"].tolist()):
        rows.append(f"1 {a!r} {l!r} {r!r}" if e else "0 0 0 0")
    case = f"output/{scenario}/{method}/{sub}"
    (ROOT / case / "error_maps.txt").write_text("\n".join(rows) + "\n")
    print(f"data/reference/{scenario}/lat.ply {case}/error_maps.txt output/mesh/{method}_{scenario}_{sub}")


def main():
    (ROOT / "output" / "mesh").mkdir(parents=True, exist_ok=True)

    for scenario in SCENARIOS_ALL_MAPS:
        sub = subs(scenario, "per_face")[-1]              # all maps
        line("per_face", scenario, sub, load(scenario, "per_face", sub))

    best = None                                           # (paired area, triplet, per-face, Roney)
    for sub in subs(SCENARIO_TRIPLET, "roney"):
        per_face, roney = load(SCENARIO_TRIPLET, "per_face", sub), load(SCENARIO_TRIPLET, "roney", sub)
        paired = roney["area"][roney["estimable"] & per_face["estimable"]].sum()
        if best is None or paired > best[0]:
            best = (paired, sub, per_face, roney)
    _, sub, per_face, roney = best
    line("per_face", SCENARIO_TRIPLET, sub, per_face)
    line("roney", SCENARIO_TRIPLET, sub, roney)


if __name__ == "__main__":
    main()
