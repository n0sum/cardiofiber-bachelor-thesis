#!/usr/bin/env python3
"""
Random per-zone tags for a scenario.

Reads data/reference/LA_um.{pts,elem,lon}, picks N random faces as centers
and assigns each face to the nearest center (Euclidean distance between centroids):
these are the N zones. Each zone gets one of the tags of CSV, at random with equal
probability. A single generator, numpy.random.default_rng(SEED), for centers and
tags. N = number of faces gives one zone per face.

Writes data/reference/NAME/LA_um: .pts and .lon copied, .elem with the same
triangles and the new tag.

Usage, from the repository root:
    python simulator/stage1/opencarp/random_tags.py NAME CSV SEED N
    python simulator/stage1/opencarp/random_tags.py aniso_zone_600 experiments/region_cv/aniso_random.csv 0 600
"""

import csv
import shutil
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

# ------------------------------------------------------------------- paths
ROOT   = Path(__file__).resolve().parents[3]                     # repository root
SRC    = ROOT / "data" / "reference" / "LA_um"                   # without extension
# --------------------------------------------------------------------------


def main():
    name, cv_csv, seed, n = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    dst = ROOT / "data" / "reference" / name / "LA_um"
    tags = [int(r["region_tag"]) for r in
            csv.DictReader(l for l in open(cv_csv) if not l.startswith("#"))]
    pts = np.loadtxt(f"{SRC}.pts", skiprows=1, dtype=np.float64)
    tri = np.loadtxt(f"{SRC}.elem", skiprows=1, usecols=(1, 2, 3), dtype=np.int64)
    bary = pts[tri].mean(1)

    rng = np.random.default_rng(seed)
    centers = rng.choice(len(tri), size=n, replace=False)
    zone = cKDTree(bary[centers]).query(bary)[1]
    tag = rng.choice(tags, size=n)[zone]

    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(f"{SRC}.pts", f"{dst}.pts")
    shutil.copyfile(f"{SRC}.lon", f"{dst}.lon")
    with open(f"{dst}.elem", "w") as f:
        f.write(f"{len(tri)}\n")
        for (a, b, c), t in zip(tri, tag):
            f.write(f"Tr {a} {b} {c} {t}\n")


if __name__ == "__main__":
    main()
