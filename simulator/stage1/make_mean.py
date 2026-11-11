import h5py, numpy as np, meshio
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]                # repository root
H5   = ROOT / "data" / "external" / "biatrial_model_completeFit.h5"
OUT  = ROOT / "data" / "stage1"

# mean shape of the bi-atrial SSM: model/mean is [x0, y0, z0, x1, ...] (matches representer/points)
with h5py.File(H5, "r") as f:
    pts  = f["model/mean"][()].reshape(-1, 3)
    tris = f["representer/cells"][()].T.astype(np.int64)
print(f"vertices: {len(pts)}   triangles: {len(tris)}")

OUT.mkdir(parents=True, exist_ok=True)
meshio.write_points_cells(OUT / "ssm_mean_surface.vtk", pts, [("triangle", tris)])
print(f"written: {OUT / 'ssm_mean_surface.vtk'}")
