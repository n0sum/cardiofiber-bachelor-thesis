#!/usr/bin/env python3
"""
Method 2 of Roney et al. 2019 on triplets of LAT maps, for the comparison with the
per-face estimator (cardiofiber/): geodesic patch, planar fit per
map, centered ellipse on the velocity vectors g / |g|^2.

With three points the conic passes exactly through all three: it is also the
solution of Roney's radial-residual fit, so the genetic algorithm is not
needed.

Usage, from the repository root:
    python estimator/roney.py experiments/roney.conf SCENARIO [SCENARIO ...]
Reads only data/reference/SCENARIO/lat.ply, never the ground truth. The
triplets are those of the lat_k maps of lat.ply; k is the site index. Writes,
next to output/SCENARIO/per_face/ of the per-face estimator,
    output/SCENARIO/roney/sub_i-j-k/estimate.ply   the estimate in format D (mesh of
                                                   lat.ply; NaN on the non-estimable faces)
    output/SCENARIO/roney/sub_i-j-k/estimate.csv   face_id,k_used
with mu = ln(lambda1 lambda2) / 4, alpha = ln(lambda1 / lambda2) / 4, lambda1 >= lambda2
the eigenvalues of M = P^-1 (CV_L^2, CV_T^2). k_used: the maps of the triplet that pass
the planar fit in the patch of the face (0-3).

Units: gradients in ms/mm, velocities in mm/ms.
"""

import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra

ROOT = Path(__file__).resolve().parents[1]                # repository root
sys.path.insert(0, str(ROOT / "simulator"))
from make_truth import write_ply                          # the writer of format D

KEYS = ["patch.radius_mm", "patch.n_points", "patch.seed", "fit.max_residual", "mask.min_singular_value"]


def read_conf(path):
    kv = {}
    for line in open(path):
        line = line.split("#")[0].strip()
        if line:
            k, v = line.split("=")
            kv[k.strip()] = v.strip()
    return {k: float(kv[k]) for k in KEYS}               # KeyError if a key is missing


def read_ply(path):
    """Binary little endian PLY with vertices and triangles: vertex properties all scalar,
    faces 'list uchar int' followed by scalar properties (double, float or int)."""
    f = open(path, "rb")
    elems, cur = {}, None
    while (line := f.readline().decode("ascii").strip()) != "end_header":
        w = line.split()
        if w[0] == "element":
            cur = w[1]
            elems[cur] = [int(w[2]), []]
        elif w[0] == "property" and w[1] == "list":
            elems[cur][1] += [("n", "u1"), ("a", "<i4"), ("b", "<i4"), ("c", "<i4")]
        elif w[0] == "property":
            elems[cur][1].append((w[2], {"double": "<f8", "float": "<f4", "int": "<i4"}[w[1]]))
    out = {}
    for name, (n, dt) in elems.items():
        out[name] = np.fromfile(f, dtype=np.dtype(dt), count=n)
    return out


def geometry(V, F, radius, n_points, seed):
    """Per face: indices of the chosen points (padded with -1 up to n_points), 2D coordinates
    (X, Y) of the points, map A (3x2) from the MDS plane to R^3, unit normal of the face."""
    nv, nf = len(V), len(F)
    e = np.concatenate([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]])
    e = np.unique(np.sort(e, axis=1), axis=0)
    w = np.linalg.norm(V[e[:, 0]] - V[e[:, 1]], axis=1)
    G = coo_matrix((np.r_[w, w], (np.r_[e[:, 0], e[:, 1]], np.r_[e[:, 1], e[:, 0]])), shape=(nv, nv)).tocsr()

    # neighbors within the radius of each vertex, in blocks so as not to keep the dense matrix
    near = [None] * nv
    for s in range(0, nv, 2000):
        idx = np.arange(s, min(s + 2000, nv))
        D = dijkstra(G, indices=idx, limit=radius)
        for r, v in enumerate(idx):
            j = np.flatnonzero(np.isfinite(D[r]))
            near[v] = (j, D[r, j])

    rng = np.random.default_rng(seed)
    pts = -np.ones((nf, n_points), dtype=np.int64)
    XY = np.zeros((nf, n_points, 2))
    A = np.zeros((nf, 3, 2))
    c_all = V[F].mean(axis=1)
    for i in range(nf):
        c = c_all[i]
        dist = np.full(nv, np.inf)
        for v in F[i]:
            j, d = near[v]
            dist[j] = np.minimum(dist[j], np.linalg.norm(c - V[v]) + d)
        patch = np.flatnonzero(dist <= radius)

        Dp = dijkstra(G[patch][:, patch])
        n = len(patch)
        J = np.eye(n) - 1.0 / n
        lam, U = np.linalg.eigh(-0.5 * J @ (Dp ** 2) @ J)
        xy = U[:, [-1, -2]] * np.sqrt(np.maximum(lam[[-1, -2]], 0))   # classical MDS, centered

        P3 = V[patch] - V[patch].mean(axis=0)
        A[i] = np.linalg.lstsq(xy, P3, rcond=None)[0].T

        sel = rng.choice(n, n_points, replace=False) if n > n_points else np.arange(n)
        pts[i, :len(sel)] = patch[sel]
        XY[i, :len(sel)] = xy[sel]

    nrm = np.cross(V[F[:, 1]] - V[F[:, 0]], V[F[:, 2]] - V[F[:, 0]])
    return pts, XY, A, nrm / np.linalg.norm(nrm, axis=1)[:, None]


def estimate(X, ok, A, nrm, min_sv):
    """X: (nf, 3, 2) velocity vectors, ok: (nf,) all three maps valid.
    Conic x^T P x = 1 through the three points; an ellipse with semi-axes CV_L, CV_T has conic M^-1, M = P^-1.
    Returns ok, the eigenvalues lambda1 >= lambda2 of M (CV_L^2, CV_T^2) and the fiber in 3D."""
    R = np.stack([X[..., 0] ** 2, np.sqrt(2) * X[..., 0] * X[..., 1], X[..., 1] ** 2], axis=-1)
    R[~ok] = np.eye(3)                                    # dummy rows: the face stays non-estimable
    ok = ok & (np.linalg.svd(R, compute_uv=False)[:, -1] > min_sv)
    R[~ok] = np.eye(3)
    p = np.linalg.solve(R, np.ones(R.shape[:2])[..., None])[..., 0]
    P = np.stack([np.stack([p[:, 0], p[:, 1] / np.sqrt(2)], -1),
                  np.stack([p[:, 1] / np.sqrt(2), p[:, 2]], -1)], -2)
    lp, U = np.linalg.eigh(P)                             # ascending; M has the same eigenvectors
    ok = ok & (lp[:, 0] > 0)
    lp = np.where(ok[:, None], lp, 1.0)
    lam, lam2, fib2 = 1 / lp[:, 0], 1 / lp[:, 1], U[:, :, 0]   # lambda1 of M = 1 / smaller eigenvalue of P
    f = np.einsum("fij,fj->fi", A, fib2)
    f -= np.sum(f * nrm, axis=1)[:, None] * nrm
    with np.errstate(divide="ignore", invalid="ignore"):
        f /= np.linalg.norm(f, axis=1)[:, None]
    return ok, lam, lam2, f


def main():
    conf = read_conf(sys.argv[1])
    radius, n_points = conf["patch.radius_mm"], int(conf["patch.n_points"])

    geo = None
    for scenario in sys.argv[2:]:
        lat = read_ply(ROOT / "data" / "reference" / scenario / "lat.ply")
        F = np.stack([lat["face"][k] for k in "abc"], axis=1).astype(np.int64)
        if geo is None:                                   # the mesh is the same in all scenarios
            V = np.stack([lat["vertex"][k] for k in "xyz"], axis=1).astype(np.float64)
            geo = geometry(V, F, radius, n_points, int(conf["patch.seed"]))
            pts, XY, A, nrm = geo
            real = pts >= 0
            D = np.concatenate([real[..., None] * 1.0, XY], axis=-1)   # rows [1, X, Y], zero if padded
            pinv = np.linalg.pinv(D)                      # (nf, 3, n_points)
            n_real = real.sum(axis=1)

        maps = [int(p.removeprefix("lat_")) for p in lat["vertex"].dtype.names if p.startswith("lat_")]
        g, valid = [], []
        for k in maps:
            t = np.where(real, lat["vertex"][f"lat_{k}"].astype(np.float64)[np.maximum(pts, 0)], 0.0)
            coef = np.einsum("fij,fj->fi", pinv, t)
            res = np.einsum("fjk,fk->fj", D, coef) - t   # zero on the padded rows
            rms = np.sqrt((res ** 2).sum(axis=1) / n_real)
            gk = coef[:, 1:]
            gn = np.linalg.norm(gk, axis=1)
            valid.append((gn > 0) & (rms <= conf["fit.max_residual"] * gn * radius))
            g.append(gk)
        g, valid = np.stack(g, axis=1), np.stack(valid, axis=1)        # (nf, maps, 2), (nf, maps)

        results = {}                                      # site triplet -> estimate
        for pos in itertools.combinations(range(len(maps)), 3):         # positions in maps
            s = list(pos)
            sub = tuple(maps[i] for i in pos)
            k_used, ok = valid[:, s].sum(axis=1), valid[:, s].all(axis=1)
            with np.errstate(divide="ignore", invalid="ignore"):
                x = g[:, s] / (g[:, s] ** 2).sum(axis=-1, keepdims=True)
            x = np.where(ok[:, None, None], x, 0.0)
            results[sub] = (k_used, *estimate(x, ok, A, nrm, conf["mask.min_singular_value"]))

        vertex = np.rec.fromarrays([lat["vertex"][k] for k in "xyz"], dtype=[("x", "<f8"), ("y", "<f8"), ("z", "<f8")])
        for sub, (k_used, ok, lam, lam2, f) in results.items():
            with np.errstate(divide="ignore", invalid="ignore"):
                mu = np.where(ok, np.log(lam * lam2) / 4, np.nan)
                alpha = np.where(ok, np.log(lam / lam2) / 4, np.nan)
            f = np.where(ok[:, None], f, np.nan)
            face = np.rec.fromarrays([np.full(len(F), 3), F[:, 0], F[:, 1], F[:, 2],
                                      f[:, 0], f[:, 1], f[:, 2], mu, alpha],
                                     dtype=[("n", "u1"), ("a", "<i4"), ("b", "<i4"), ("c", "<i4"),
                                            ("fx", "<f8"), ("fy", "<f8"), ("fz", "<f8"),
                                            ("mu", "<f8"), ("alpha", "<f8")])
            out = ROOT / "output" / scenario / "roney" / ("sub_" + "-".join(map(str, sub)))
            out.mkdir(parents=True, exist_ok=True)
            write_ply(out / "estimate.ply", vertex, face)
            pd.DataFrame({"face_id": np.arange(len(F)), "k_used": k_used}).to_csv(out / "estimate.csv", index=False)


if __name__ == "__main__":
    main()
