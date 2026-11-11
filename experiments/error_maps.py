#!/usr/bin/env python3
"""
Error maps in PLY for MeshLab, from the errors.csv files of the grid.
For each case writes next to its errors.csv the file error_maps.txt, only
numbers: one row per face, in face_id order, "estimable angle cv_l ratio",
0 0 0 0 on the non-estimable faces. Then prints a line "MESH NUMBERS PREFIX",
to be passed to ErrorMaps, which writes PREFIX_{angle,cv_l,ratio}.ply.
MESH is the LAT PLY of the scenario, NUMBERS the error_maps.txt, PREFIX
output/mesh/METHOD_SCENARIO_SUB. The cases:
    per-face, all maps, on SCENARIO_ALL_MAPS and on SCENARIO_TRIPLET
    per-face and Roney on the triplet of SCENARIO_TRIPLET with the largest
    paired area: faces estimable by both methods, area-weighted

Usage, from the repository root:
    python experiments/error_maps.py | while read -r mesh numbers prefix; do
        ./estimator/build/ErrorMaps "$mesh" "$numbers" "$prefix"
    done
"""

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]                # repository root
sys.path.insert(0, str(ROOT / "analysis"))
from analysis import load                                 # the error columns of errors.csv

SCENARIO_ALL_MAPS = "aniso_zone_1"
SCENARIO_TRIPLET = "aniso_zone_600"


def line(method, scenario, sub):
    """Writes error_maps.txt of the case and prints its line. Numbers with repr, which reads back
    exactly as a double; no "nan", which >> in C++ does not read."""
    case = ROOT / "output" / scenario / method / sub
    df = load(case / "errors.csv")
    rows = []
    for e, a, l, r in zip(df["estimable"].tolist(), df["angle"].tolist(), df["cv_l"].tolist(),
                          df["ratio"].tolist()):
        rows.append(f"1 {a!r} {l!r} {r!r}" if e == 1 else "0 0 0 0")
    (case / "error_maps.txt").write_text("\n".join(rows) + "\n")
    numbers = f"output/{scenario}/{method}/{sub}/error_maps.txt"
    print(f"data/reference/LA_LAT_{scenario}.ply {numbers} output/mesh/{method}_{scenario}_{sub}")


def main():
    (ROOT / "output" / "mesh").mkdir(parents=True, exist_ok=True)

    for scenario in (SCENARIO_ALL_MAPS, SCENARIO_TRIPLET):
        subs = [d.name for d in (ROOT / "output" / scenario / "per_face").glob("sub_*")]
        line("per_face", scenario, max(subs, key=lambda s: s.count("-")))   # all maps

    paired = {}                                           # triplet -> paired area
    for d in (ROOT / "output" / SCENARIO_TRIPLET / "roney").glob("sub_*"):
        roney = pd.read_csv(d / "errors.csv")
        per_face = pd.read_csv(ROOT / "output" / SCENARIO_TRIPLET / "per_face" / d.name / "errors.csv")
        both = (roney["estimable"] == 1) & (per_face["estimable"] == 1)
        paired[d.name] = roney["area"][both].sum() / roney["area"].sum()
    sub = max(paired, key=paired.get)
    line("per_face", SCENARIO_TRIPLET, sub)
    line("roney", SCENARIO_TRIPLET, sub)


if __name__ == "__main__":
    main()
