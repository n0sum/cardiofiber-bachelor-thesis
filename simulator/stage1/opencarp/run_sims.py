#!/usr/bin/env python3
"""
One openCARP eikonal simulation for each row of pacing_sites.csv.

The stimulus point is the centroid of cell cell_id on the mesh in um,
rounded to the nearest integer. The sites are defined by cell
index, not by coordinates. Row k of the CSV is map k.

The CVs come from the scenario CSV (region_tag,CV_L,CV_T in mm/ms), the single
source: here they are converted to mm/s for openCARP, one gregion per tag. The
same tags give the single eikonal domain (phys_region).

MESH_BASE is the mesh in um without extension: its per-face tags must
be those of the CSV.

The runs end up in data/reference/sims/NAME/run_k/. openCARP has no
option to overwrite an existing run, so the folder is
deleted first.

Usage, from the repository root:
    python simulator/stage1/opencarp/run_sims.py NAME CSV MESH_BASE
    python simulator/stage1/opencarp/run_sims.py aniso_zone_70 experiments/region_cv/aniso_random.csv data/reference/aniso_zone_70/LA_um
"""

import csv
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

# ------------------------------------------------------------------- paths
HERE     = Path(__file__).resolve().parent
ROOT     = HERE.parents[2]                                    # repository root
SITES    = HERE.parent / "pacing_sites.csv"
SIMS_DIR = ROOT / "data" / "reference" / "sims"
# --------------------------------------------------------------------------


def region_args(cv_csv):
    """openCARP phys_region and gregion arguments from the CV CSV. The only
    place of the mm/ms -> mm/s conversion."""
    rows = list(csv.DictReader(l for l in open(cv_csv) if not l.startswith("#")))
    args = ["-num_phys_regions", "1",
            "-phys_region[0].name", "Eikonal domain",
            "-phys_region[0].ptype", "2",
            "-phys_region[0].num_IDs", str(len(rows))]
    for i, row in enumerate(rows):
        args += [f"-phys_region[0].ID[{i}]", row["region_tag"]]
    args += ["-num_gregions", str(len(rows))]
    for i, row in enumerate(rows):
        vel_l = float(row["CV_L"]) * 1000                     # mm/ms -> mm/s
        vel_t = float(row["CV_T"]) * 1000
        args += [f"-gregion[{i}].num_IDs", "1",
                 f"-gregion[{i}].ID[0]", row["region_tag"],
                 f"-gregion[{i}].dream.vel_l", f"{vel_l:g}",
                 f"-gregion[{i}].dream.vel_t", f"{vel_t:g}",
                 f"-gregion[{i}].dream.vel_n", f"{vel_t:g}"]
    return args


def main():
    name = sys.argv[1]
    cv_csv = Path(sys.argv[2]).resolve()
    mesh_um = Path(sys.argv[3]).resolve()
    sims_dir = SIMS_DIR / name
    regions = region_args(cv_csv)
    pts = np.loadtxt(f"{mesh_um}.pts", skiprows=1, dtype=np.float64)
    tri = np.loadtxt(f"{mesh_um}.elem", skiprows=1, usecols=(1, 2, 3), dtype=np.int64)
    sims_dir.mkdir(parents=True, exist_ok=True)

    for k, row in enumerate(csv.DictReader(open(SITES))):
        p = np.round(pts[tri[int(row["cell_id"])]].mean(0)).astype(np.int64)
        shutil.rmtree(sims_dir / f"run_{k}", ignore_errors=True)
        subprocess.run([HERE / "run_eik.sh", str(mesh_um), str(k), *map(str, p), *regions],
                       cwd=sims_dir, check=True)


if __name__ == "__main__":
    main()
