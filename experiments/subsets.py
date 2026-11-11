#!/usr/bin/env python3
"""
Generates one .conf for each subset of k >= 3 LAT maps of a scenario
(with six maps 20 + 15 + 6 + 1 = 42). The number of maps is the number of
rows of pacing_sites.csv (row k is map k).

Each subset is named sub_ + indices joined by hyphens (e.g. sub_0-2-4): the name of
its .conf and of its output folder.

Usage, from the repository root:
    python experiments/subsets.py SCENARIO
The .conf files go to output/SCENARIO/confs/, the outputs of Estimate to
output/SCENARIO/per_face/sub_.../errors.csv; the folder of each errors.csv is
created here.
"""

import csv
import itertools
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]                # repository root
MODEL = ROOT / "experiments" / "template.conf"
SITES = ROOT / "simulator" / "stage1" / "pacing_sites.csv"
MIN_MAPS = 3                                              # kMinGradients in cardiofiber/src/Mesh.cpp


def write_conf(path, scenario, maps, output_csv):
    """Writes a .conf to path: the thresholds of template.conf plus scenario, maps and output.
    Also creates the folder of output_csv: Estimate does not create it."""
    keys = {"lat_ply": f"data/reference/LA_LAT_{scenario}.ply",
            "ground_truth_ply": f"data/reference/LA_GT_{scenario}.ply",
            "maps": maps,
            "output_csv": output_csv}
    lines = [l for l in MODEL.read_text().splitlines() if l.strip() and not l.startswith("#")]
    lines += [f"{k} = {v}" for k, v in keys.items()]
    path.parent.mkdir(parents=True, exist_ok=True)
    (ROOT / output_csv).parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def main():
    scenario = sys.argv[1]
    out_dir = f"output/{scenario}"
    maps = range(len(list(csv.DictReader(open(SITES)))))
    for k in range(MIN_MAPS, len(maps) + 1):
        for subset in itertools.combinations(maps, k):
            sub_name = "sub_" + "-".join(map(str, subset))
            write_conf(ROOT / out_dir / "confs" / f"{sub_name}.conf", scenario,
                       ",".join(map(str, subset)), f"{out_dir}/per_face/{sub_name}/errors.csv")


if __name__ == "__main__":
    main()
