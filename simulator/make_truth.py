#!/usr/bin/env python3
"""
Ground-truth conduction field of a scenario, in format D.

    data/reference/SCENARIO/truth.ply   per vertex x y z (mm); per face the
                                        triangle, the unit fiber fx fy fz and
                                        mu, alpha (log form of CV_L, CV_T)

mu = (ln CV_L + ln CV_T) / 2, alpha = (ln CV_L - ln CV_T) / 2, CV in mm/ms,
natural log. Mesh and fibers come from the AugmentA output in data/augmenta/,
read as text: pyvista or meshio may reorder the cells and misalign faces and
fibers without errors. The faces keep the order of the .elem but are flipped,
(a, b, c) -> (a, c, b), so that the normals point outwards: this is the only
place in the pipeline where faces are flipped.

Scenarios:
    uniform   CV_L = 1.05, CV_T = 0.525 mm/ms on every face

Usage, from the repository root:
    python simulator/make_truth.py uniform
"""

import sys
from pathlib import Path

import numpy as np

ROOT     = Path(__file__).resolve().parents[1]                    # repository root
AUGMENTA = ROOT / "data" / "augmenta" / "LA_endo_with_fiber"      # without extension

COMMENT = "mu = (ln CV_L + ln CV_T) / 2, alpha = (ln CV_L - ln CV_T) / 2, CV in mm/ms, natural log"


def write_ply(path, vertex, face):
    """vertex and face are structured arrays; face starts with n, a, b, c, the list of vertices."""
    ply_type = {np.dtype("<f4"): "float", np.dtype("<f8"): "double", np.dtype("<i4"): "int"}
    header = ["ply", "format binary_little_endian 1.0", f"comment {COMMENT}", f"element vertex {len(vertex)}"]
    header += [f"property {ply_type[vertex.dtype[p]]} {p}" for p in vertex.dtype.names]
    header += [f"element face {len(face)}", "property list uchar int vertex_indices"]
    header += [f"property {ply_type[face.dtype[p]]} {p}" for p in face.dtype.names[4:]]
    header += ["end_header"]
    with open(path, "wb") as f:
        f.write(("\n".join(header) + "\n").encode("ascii"))
        vertex.tofile(f)
        face.tofile(f)


def main():
    if len(sys.argv) != 2:
        raise SystemExit("usage: python simulator/make_truth.py SCENARIO")
    scenario = sys.argv[1]
    if scenario != "uniform":
        raise SystemExit(f"unknown scenario '{scenario}'; available: uniform")

    pts = np.loadtxt(f"{AUGMENTA}.pts", skiprows=1)
    tri = np.loadtxt(f"{AUGMENTA}.elem", skiprows=1, usecols=(1, 2, 3), dtype=np.int64)
    fib = np.loadtxt(f"{AUGMENTA}.lon", skiprows=1)[:, :3]          # fiber; the sheet is not needed

    tri = tri[:, [0, 2, 1]]                                          # AugmentA normals point inwards
    fib = fib / np.linalg.norm(fib, axis=1)[:, None]

    cv_l = np.full(len(tri), 1.05)                                   # mm/ms
    cv_t = np.full(len(tri), 0.525)
    mu = (np.log(cv_l) + np.log(cv_t)) / 2
    alpha = (np.log(cv_l) - np.log(cv_t)) / 2

    vertex = np.rec.fromarrays([pts[:, 0], pts[:, 1], pts[:, 2]],
                               dtype=[("x", "<f8"), ("y", "<f8"), ("z", "<f8")])
    face = np.rec.fromarrays([np.full(len(tri), 3), tri[:, 0], tri[:, 1], tri[:, 2],
                              fib[:, 0], fib[:, 1], fib[:, 2], mu, alpha],
                             dtype=[("n", "u1"), ("a", "<i4"), ("b", "<i4"), ("c", "<i4"),
                                    ("fx", "<f8"), ("fy", "<f8"), ("fz", "<f8"),
                                    ("mu", "<f8"), ("alpha", "<f8")])

    out = ROOT / "data" / "reference" / scenario
    out.mkdir(parents=True, exist_ok=True)
    write_ply(out / "truth.ply", vertex, face)


if __name__ == "__main__":
    main()
