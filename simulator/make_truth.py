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
    painted   alpha and mu painted by hand in MeshLab, in grayscale, on the
              AugmentA mesh with flipped faces: simulator/painted_alpha.ply and
              simulator/painted_mu.ply. Only their per-vertex colors are read,
              matched by vertex index; g = mean of R, G, B / 255 (0 black,
              1 white), the alpha channel is ignored.
                alpha = 1/2 ln 1.6 + g (1/2 ln 2.4 - 1/2 ln 1.6)   r from 1.6 to 2.4
                mu    = -0.6 + 0.6 g                               e^mu from 0.55 to 1.0 mm/ms
              Per face the mean of the values of its three vertices.
    r1.2, r1.41, r1.6, r2.4, r3
              the sweep in r: mu of uniform, constant, and alpha = 1/2 ln r
              constant, with r = 1.2, sqrt(2) exactly, 1.6, 2.4, 3; r = 2 is
              uniform itself. Same fibers and mu as uniform: only the
              flattening of the ellipse changes. Above r = sqrt(2) the pedal
              curve of the front ellipse (the curve of the velocity vectors)
              is no longer convex.

Usage, from the repository root:
    python simulator/make_truth.py SCENARIO     # e.g. uniform, painted, r1.41
"""

import sys
from pathlib import Path

import meshio
import numpy as np

ROOT     = Path(__file__).resolve().parents[1]                    # repository root
AUGMENTA = ROOT / "data" / "augmenta" / "LA_endo_with_fiber"      # without extension
PAINTED  = ROOT / "simulator"                                     # painted_alpha.ply, painted_mu.ply
MU_UNIFORM = (np.log(1.05) + np.log(0.525)) / 2                    # mu of uniform, also of the sweep in r
SWEEP_ALPHA = {"r1.2": np.log(1.2) / 2, "r1.41": np.log(2) / 4,     # alpha = 1/2 ln r; r1.41 is r = sqrt(2)
               "r1.6": np.log(1.6) / 2, "r2.4": np.log(2.4) / 2, "r3": np.log(3) / 2}
SCENARIOS = ("uniform", "painted", *SWEEP_ALPHA)

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


def gray(path):
    """Gray level per vertex of a PLY painted in MeshLab: mean of R, G, B / 255, 0 black, 1 white.
    meshio reads the uchar colors of a binary PLY as int8 (255 -> -1): they are reinterpreted as
    uint8, which leaves those of a text PLY unchanged."""
    data = meshio.read(path).point_data
    return np.mean([data[c].astype(np.uint8) for c in ("red", "green", "blue")], axis=0) / 255


def painted(tri):
    """mu and alpha per face of the painted scenario: per vertex from the gray levels, then the mean of
    the three vertices of each face."""
    g_alpha, g_mu = gray(PAINTED / "painted_alpha.ply"), gray(PAINTED / "painted_mu.ply")
    alpha_lo, alpha_hi = np.log(1.6) / 2, np.log(2.4) / 2          # black r = 1.6, white r = 2.4
    alpha = alpha_lo + g_alpha * (alpha_hi - alpha_lo)
    mu = -0.6 + 0.6 * g_mu                                           # black e^mu = 0.55, white 1.0 mm/ms
    return mu[tri].mean(axis=1), alpha[tri].mean(axis=1)


def main():
    if len(sys.argv) != 2:
        raise SystemExit("usage: python simulator/make_truth.py SCENARIO")
    scenario = sys.argv[1]
    if scenario not in SCENARIOS:
        raise SystemExit(f"unknown scenario '{scenario}'; available: {', '.join(SCENARIOS)}")

    pts = np.loadtxt(f"{AUGMENTA}.pts", skiprows=1)
    tri = np.loadtxt(f"{AUGMENTA}.elem", skiprows=1, usecols=(1, 2, 3), dtype=np.int64)
    fib = np.loadtxt(f"{AUGMENTA}.lon", skiprows=1)[:, :3]          # fiber; the sheet is not needed

    tri = tri[:, [0, 2, 1]]                                          # AugmentA normals point inwards
    fib = fib / np.linalg.norm(fib, axis=1)[:, None]

    if scenario == "uniform":
        cv_l = np.full(len(tri), 1.05)                               # mm/ms
        cv_t = np.full(len(tri), 0.525)
        mu = (np.log(cv_l) + np.log(cv_t)) / 2
        alpha = (np.log(cv_l) - np.log(cv_t)) / 2
    elif scenario == "painted":
        mu, alpha = painted(tri)
    else:
        mu = np.full(len(tri), MU_UNIFORM)
        alpha = np.full(len(tri), SWEEP_ALPHA[scenario])

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
