#!/usr/bin/env python3
"""
Error maps in PLY for MeshLab, from the runs of the grid. The errors come from
load() of analysis/analysis.py, against the ground truth.
For each case writes next to its estimate.ply the file error_maps.txt, only
numbers: one row per face, in face_id order, "estimable angle mean_cv ratio",
0 0 0 0 on the non-estimable faces. Then creates the folder of the case,
CASEDIR = output/mesh/METHOD_SCENARIO_SUB, and prints a line
"MESH NUMBERS CASEDIR": MESH is data/reference/SCENARIO/lat.ply, NUMBERS the
error_maps.txt. run_all.sh passes them to ErrorMaps with PREFIX CASEDIR/error,
which writes CASEDIR/error_{angle,mean_cv,ratio}.ply. The cases:
    per-face, all maps, on each scenario of SCENARIOS_ALL_MAPS
    per-face and Roney on the triplet of SCENARIO_TRIPLET with the largest
    paired area: faces estimable by both methods, area-weighted

With NOISE (e.g. 0.1, NOISES in run_all.sh) the same cases on the scenarios
SCENARIO_noiseNOISE: SCENARIO in the paths and in CASEDIR has the suffix, e.g.
output/mesh/per_face_painted_noise0.1_sub_0-1-2-3-4-5/.

Usage, from the repository root:
    python experiments/error_maps.py [NOISE] | while read -r mesh numbers case_dir; do
        ./estimator/build/ErrorMaps "$mesh" "$numbers" "$case_dir/error"
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
    """Writes error_maps.txt of the case, creates its folder in output/mesh/ and prints its line. Numbers with repr, which reads back
    exactly as a double; no "nan", which >> in C++ does not read."""
    rows = []
    for e, a, m, r in zip(df["estimable"].tolist(), df["angle"].tolist(), df["mean_cv"].tolist(),
                          df["ratio"].tolist()):
        rows.append(f"1 {a!r} {m!r} {r!r}" if e else "0 0 0 0")
    case = f"output/{scenario}/{method}/{sub}"
    (ROOT / case / "error_maps.txt").write_text("\n".join(rows) + "\n")
    case_dir = f"output/mesh/{method}_{scenario}_{sub}"
    (ROOT / case_dir).mkdir(parents=True, exist_ok=True)
    print(f"data/reference/{scenario}/lat.ply {case}/error_maps.txt {case_dir}")


def main():
    suffix = f"_noise{sys.argv[1]}" if len(sys.argv) > 1 else ""   # NOISE: the scenarios SCENARIO_noiseNOISE
    for scenario in SCENARIOS_ALL_MAPS:
        scenario += suffix
        sub = subs(scenario, "per_face")[-1]              # all maps
        line("per_face", scenario, sub, load(scenario, "per_face", sub))

    triplet_scenario = SCENARIO_TRIPLET + suffix
    best = None                                           # (paired area, triplet, per-face, Roney)
    for sub in subs(triplet_scenario, "roney"):
        per_face, roney = load(triplet_scenario, "per_face", sub), load(triplet_scenario, "roney", sub)
        paired = roney["area"][roney["estimable"] & per_face["estimable"]].sum()
        if best is None or paired > best[0]:
            best = (paired, sub, per_face, roney)
    _, sub, per_face, roney = best
    line("per_face", triplet_scenario, sub, per_face)
    line("roney", triplet_scenario, sub, roney)


if __name__ == "__main__":
    main()
