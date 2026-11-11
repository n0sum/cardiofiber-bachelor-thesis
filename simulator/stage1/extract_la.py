import pyvista as pv, numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]                # repository root
DIR  = ROOT / "data" / "stage1"
MASK = ROOT / "simulator" / "stage1" / "la_vertex_mask.npy"

# Triangle orientation: AugmentA flips the orientation of all triangles with respect
# to the SSM (elem == ssm_tri[:, ::-1]), same set and same order. The
# reference connectivity for all shapes is LA_endo_with_fiber.elem,
# never representer/cells.

LA_REGION = 0   # RegionId of the LA among the two connected components (atria) of the surface

mesh = pv.read(DIR / "ssm_mean_surface.vtk")

# stamp each vertex with its original index (0..39868) BEFORE splitting
mesh.point_data["orig_id"] = np.arange(mesh.n_points, dtype=np.int64)

conn = mesh.connectivity()
rid  = np.asarray(conn.cell_data["RegionId"])

la_cells = np.where(rid == LA_REGION)[0]
la_surf  = conn.extract_cells(la_cells).extract_surface(algorithm="dataset_surface")   # pyvista 0.48.4 default, pinned

# the original indices survive in the field we created ourselves
la_point_ids = np.asarray(la_surf.point_data["orig_id"]).astype(np.int64)
np.save(MASK, np.sort(la_point_ids))

la_surf.save(DIR / "ssm_mean_LA.vtk")
print(f"LA: {la_surf.n_points} vertices, {la_surf.n_cells} triangles")
print(f"mask saved: la_vertex_mask.npy ({la_point_ids.size} indices)")
print(f"index range: [{la_point_ids.min()}, {la_point_ids.max()}]  expected max < 39869")
