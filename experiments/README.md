# Experiments

The pipeline, for each scenario (`ALL_SCENARIOS` in `run_all.sh`, built from
`NOISES`, `RATIOS` and `CORE_NOISE` at its top; in this order `painted`,
`painted_noise0.05`, `painted_noise0.1`, `painted_noise0.3`,
`painted_r1.2_noise0.05`, `painted_r2_noise0.05`, `painted_r2.8_noise0.05`):

```
simulator/make_truth.py TRUTH            ->  data/reference/TRUTH/truth.ply      ground truth, format D
Simulate truth.ply SITES lat.ply SIGMA   ->  data/reference/SCENARIO/lat.ply     six LAT maps
experiments/subsets.py SCENARIO          ->  output/SCENARIO/confs/sub_….conf
Estimate output/SCENARIO/confs/*         ->  output/SCENARIO/per_face/sub_…/estimate.{ply,csv}
estimator/roney.py roney.conf SC..       ->  output/SCENARIO/roney/sub_i-j-k/estimate.{ply,csv}
analysis/analysis.py                     ->  output/figures/
experiments/error_maps.py + ErrorMaps + FiberView -> output/mesh/
```

## Noise on the LATs

`Simulate` takes a fourth argument SIGMA, in ms. After computing each map it
adds to every vertex (also those of the stimulus sphere) a number drawn from a
Gaussian with mean 0 and standard deviation SIGMA: additive noise in ms,
independent per vertex and per map, from one `std::mt19937` with a fixed seed
(`NOISE_SEED` in `Simulate.cpp`). With SIGMA = 0 nothing is drawn. Since the
seed is the same for every call, the noise at σ = 0.3 is exactly 6 times that
at σ = 0.05, and the scenarios of the core sample have the same realization of
the noise as `painted_noise0.05`.

`run_all.sh` writes the ground truths `painted`, `painted_r1.2`, `painted_r2`
and `painted_r2.8`, simulates `painted` clean (SIGMA = 0) and then the noisy
scenarios, each a folder with a copy of a `truth.ply` and the LATs with SIGMA:
`painted_noiseN` (ground truth of `painted`) for each N of `NOISES`
(`0.05 0.1 0.3`), and `painted_rR_noise0.05` (ground truth of `painted_rR`,
SIGMA = `CORE_NOISE`) for each R of `RATIOS` (`1.2 2 2.8`). The folders
`data/reference/painted_rR/` hold only the ground truth. `analysis.py` and
`error_maps.py` run once, without arguments.

The two estimators read only `lat.ply` and write their estimate in format D:
they never see the ground truth. The errors are computed only in the analysis,
by `load()` of `analysis/analysis.py`, with the same formulas for the two
methods.

Units: lengths in mm, times in ms, velocities in mm/ms.

## Scenarios

`simulator/make_truth.py TRUTH`; all use the AugmentA mesh and fibers and
differ only in μ and α:

- `painted`: μ and α painted by hand in MeshLab, in grayscale, in
  `simulator/painted_mu.ply` and `simulator/painted_alpha.ply` (the AugmentA
  mesh with flipped faces). Only the per-vertex colors are read, matched by
  vertex index. Per vertex g = (R + G + B) / 3 / 255, 0 black and 1 white (the
  alpha channel is ignored); then
  α = ½ ln 1.6 + g (½ ln 2.4 − ½ ln 1.6), r = e^(2α) from 1.6 (black) to
  2.4 (white), and μ = −0.6 + 0.6 g, e^μ from 0.55 (black) to 1.0 mm/ms
  (white). Per face, the mean of the values (not of the colors) of its three
  vertices.
- `painted_r1.2`, `painted_r2`, `painted_r2.8`: μ of `painted` (the same
  function, identical face by face) and α = ½ ln r constant on every face,
  r = 1.2, 2, 2.8. r is constant on purpose: a varying r would mix with the
  position on the mesh. Only the flattening of the ellipse changes.

The roles of the seven scenarios in the analysis:

| role | scenario | ground truth | σ (ms) |
|---|---|---|---|
| grid | `painted` | `painted` | 0 |
| grid | `painted_noise0.05` | `painted` | 0.05 |
| grid | `painted_noise0.1` | `painted` | 0.1 |
| grid | `painted_noise0.3` | `painted` | 0.3 |
| core sample | `painted_r1.2_noise0.05` | `painted_r1.2` | 0.05 |
| core sample | `painted_r2_noise0.05` | `painted_r2` | 0.05 |
| core sample | `painted_r2.8_noise0.05` | `painted_r2.8` | 0.05 |

## Format D: conduction field

The ground truth and the estimates of both methods. Text PLY, one line per
vertex and one per face (`3 a b c fx fy fz mu alpha`), doubles with 17
significant digits, which read back exactly:

```
ply
format ascii 1.0
comment mu = (ln CV_L + ln CV_T) / 2, alpha = (ln CV_L - ln CV_T) / 2, CV in mm/ms, natural log
element vertex 21811
property double x
property double y
property double z
element face 43270
property list uchar int vertex_indices
property double fx
property double fy
property double fz
property double mu
property double alpha
end_header
```

- Back from the log form: CV_L = e^(μ+α), CV_T = e^(μ−α), r = CV_L/CV_T = e^(2α).
- Vertices in mm. Faces with outward normals: `make_truth.py` flips the
  AugmentA triangles once, `(a, b, c) -> (a, c, b)`; after that nothing flips
  them again, and every other file has the faces of `truth.ply`.
- Fiber `fx fy fz`: unit, tangent to the face; its sign is irrelevant.
- Faces without a value (non-estimable): NaN in `fx fy fz mu alpha`, written
  `nan` (`-nan` by happly when the sign bit is set). Never zeros: μ = 0 means
  1 mm/ms. happly reads them with `strtod` (see `estimator/THIRD_PARTY.md`).
- The files written by `Estimate` (happly) have one more comment line,
  `Written with hapPLY`.

## `lat.ply`: activation maps

Text PLY (`format ascii 1.0`) written by `Simulate`: per vertex `x y z lat_0 … lat_5` (double;
mm, ms), faces `vertex_indices` as in `truth.ply`. Map k is paced from the face
in row k of `simulator/stage1/pacing_sites.csv`. In the noisy scenarios
the LATs carry the Gaussian noise described above.

## Estimator `.conf`

One `.conf` file per experiment of the per-face estimator. `Estimate` takes one
or more and processes them in sequence:

```bash
./estimator/build/Estimate output/painted/confs/sub_0-1-2.conf output/painted/confs/sub_0-1-3.conf
```

`key = value`, one per line, with spaces around `=` and no spaces in the
value. Empty lines and lines starting with `#` are skipped. All 5 keys are
required, without defaults:

| key                         | meaning                                         |
|-----------------------------|-------------------------------------------------|
| `lat_ply`                   | mesh + LAT maps: the only input of the estimator |
| `maps`                      | comma-separated map indices, no spaces: `0,2,4` |
| `mask.min_singular_value`   | threshold on the smallest singular value of the per-face system |
| `output_ply`                | where to write the estimate, format D            |
| `output_csv`                | where to write `face_id,cond`                    |

The `.conf` files are generated by `subsets.py SCENARIO`, with `write_conf`:
the lines of `template.conf` (the threshold, the same for all) plus `lat_ply`,
`maps`, `output_ply` and `output_csv`. They are not in git: the recipe is the
script. `subsets.py` writes one `.conf` for each subset of at least 3 maps
(the rows of `pacing_sites.csv`; with six maps 20 + 15 + 6 + 1 = 42) in
`output/SCENARIO/confs/`, and creates the output folder
`output/SCENARIO/per_face/ID/`. The ID is `sub_` with the site indices joined
by hyphens, e.g. `sub_0-2-4`, also the name of the `.conf`. The six maps are
`sub_0-1-2-3-4-5`. A `.conf` written by hand works too, but the folder of the
outputs must already exist: `Estimate` does not create it.

## Non-estimable faces (per-face estimator)

Each face uses all the maps of `maps`, one row per map in its linear system.
It is non-estimable if it has fewer than 3 maps, if the smallest singular value of the
system does not exceed `mask.min_singular_value` (nearly parallel fronts) or if
the smaller eigenvalue of the estimated tensor (CV_T^2) does not exceed 1e-8 times
the larger one (CV_L^2): a test of numerical zero, CV_T / CV_L > 1e-4, not a physical
threshold (`kMinEigenRatio` in `cardiofiber/src/Mesh.cpp`). In the last case
`cond` is computed but the estimate is NaN.

## `estimate.csv`

One row per face, the quantities of the estimator that are not in format D:

| method    | columns          | meaning |
|-----------|------------------|---------|
| per-face  | `face_id,cond`   | `cond` = `sigma_max / sigma_min` of the system of the face; empty where it is not computed (fewer than 3 maps, or smallest singular value below threshold) |
| Roney     | `face_id,k_used` | maps of the triplet that pass the planar fit in the patch of the face (0-3). With `k_used < 3` the face is not estimable because of the fit; with `k_used = 3` and NaN in the estimate the conic is not an ellipse or the system is nearly singular |

Scenario, method and subset are given by the folder,
`output/SCENARIO/METHOD/sub_…/`, not by the columns.

## Roney's method

`python estimator/roney.py experiments/roney.conf SCENARIO ...` applies method 2
of Roney et al. 2019 to all the 20 map triplets of each scenario and writes
`output/SCENARIO/roney/sub_i-j-k/estimate.{ply,csv}`. The keys of `roney.conf`
are documented in the file, the method and its departures from the paper in
`estimator/README_roney.md`.

## The experiment grid

```
output/
  SCENARIO/                           one per scenario
    confs/sub_….conf                  subsets.py, 42 subsets
    per_face/sub_…/estimate.ply       Estimate, 42 subsets
    per_face/sub_…/estimate.csv
    roney/sub_i-j-k/estimate.ply      estimator/roney.py, 20 triplets
    roney/sub_i-j-k/estimate.csv
  figures/                            analysis/analysis.py, tables and PNGs
  mesh/                               all the PLYs for MeshLab (ErrorMaps, FiberView), below
```

In the folders of the error-map cases (`error_maps.py`) there is also
`error_maps.txt`. `output/mesh/`, for each scenario S of the grid (`painted`,
`painted_noise0.05`, `painted_noise0.1`, `painted_noise0.3`), with TRIPLET the
triplet chosen by `error_maps.py` (below), the same for every S:

```
output/mesh/
  per_face_S_sub_0-1-2-3-4-5/   error.ply  estimate_ratio.ply  estimate_fibers.ply
  per_face_S_TRIPLET/           (the same three)
  roney_S_TRIPLET/              (the same three)
  truth_S/                      truth_ratio.ply  truth_fibers.ply
```

`run_all.sh` deletes `output/` before recomputing it, so no old results end up
in the figures.

## Analysis

`python analysis/analysis.py`, without arguments, writes everything to
`output/figures/`. Constants at its top: `GRID`, the scenarios of the grid with
their σ in ms (`painted` 0, `painted_noise0.05`, `painted_noise0.1`,
`painted_noise0.3`); `CORE`, the scenarios of the core sample with their r
(1.2, 2, 2.8); `CORE_K = 4` and `CORE_NOISE = 0.05`, its single condition; and
`SCENARIOS`, those of `GRID` and then those of `CORE`, in the order of
`run_all.sh`, for `per_face.csv` and `comparison.csv`.

`load(scenario, method, sub)` reads `data/reference/SCENARIO/truth.ply` (once per
scenario), `estimate.ply` and `estimate.csv` of the run, and returns per face:

| column      | meaning |
|-------------|---------|
| `area`      | area of the face from the vertices of the ground truth, `|(b − a) × (c − a)| / 2`, mm² |
| `estimable` | μ of the estimate finite |
| `cond` / `k_used` | from `estimate.csv` |
| `error`     | log-Euclidean distance between estimated and true tensor, below |

The tensor of a face is M = CV_L² f fᵀ + CV_T² t tᵀ (f fiber, t transverse
direction in the plane of the face). Its matrix logarithm has the same
directions and the logarithms of the eigenvalues, log M = 2μ I + 2α (f fᵀ − t tᵀ):
the three numbers of format D. The error (Arsigny et al. 2006) is

```
error = ‖log M̂ − log M‖_F
      = 2√2 · sqrt((μ̂ − μ)² + α̂² + α² − 2 α̂ α cos 2θ)
```

with θ the angle between estimated and true fiber, computed without matrix
logarithms from format D: c = |f̂·f|, cos 2θ = 2c² − 1. It is dimensionless;
with exact fiber and anisotropy, a relative error ε on both speeds gives about
2√2 ε (0.05 ≈ 1.8% on the speeds). With α = 0
(isotropic tissue) the fiber does not enter it; it does not depend on the unit
of the speeds (μ shifts by a constant); it is the same for speeds or
slownesses (log M⁻¹ = −log M).

The error is NaN on the non-estimable faces. Statistics are area-weighted:
q50, q90 and q99 (`error_q50/q90/q99` in the tables) are the values below
which 50, 90 and 99% of the area lies.

- `per_face.csv`: one row per scenario and subset of the per-face
  estimator: non-estimable area, errors on the estimable faces of the
  subset, `cond_q50/q90/q99` on the faces with computed cond.
- `comparison.csv`: one row per scenario, triplet and method (`per_face`,
  `roney`), with the non-estimable area of the method alone, the paired area
  (faces estimable by both methods with that triplet) and the errors on the
  paired faces. For the Roney rows the non-estimable area is also split in
  two: `fit_area_%` (`k_used < 3`, lost to the planar fit) and `conic_area_%`
  (`k_used = 3` and estimate NaN: the conic is not an ellipse or the system is
  nearly singular); both empty on the per-face rows. The rows with `sub` =
  `pooled` gather the paired faces of the 20 triplets of the scenario; their
  areas are means over the triplets.
- `grid.csv`: one row per (σ, k) of the scenarios of `GRID`, columns
  `scenario, noise_ms, k, subsets, non_estimable_area_%, error_mean,
  error_q50, error_q90, error_q99`.
- `core.csv`: one row per scenario of `CORE`, at k = `CORE_K`, the same
  columns with `r` in place of `noise_ms`.

A cell (scenario, k) of `grid.csv` and `core.csv` (`cell()`) takes all the
subsets with k maps of the scenario (20 triplets, 15 quadruplets, 6
quintuplets, 1 sextuplet) and concatenates their DataFrames from `load()`:
each face appears once per subset. `error_mean` = Σ area·E / Σ area and
`error_q50/q90/q99` (area-weighted quantiles) are on the estimable rows;
`non_estimable_area_%` = area of the non-estimable rows / area of all the
rows × 100, i.e. the mean over the subsets of the non-estimable area of each;
`subsets` = number of subsets of the cell. The cond (`cond_q50/q90/q99` of
`per_face.csv`, from `estimate.csv`) is in no figure.

| file | content |
|------|---------|
| `plot_grid.png` | the grid: a 4 × 4 colored table, k = 3…6 on the horizontal axis, σ = 0, 0.05, 0.1, 0.3 ms on the vertical axis as equispaced categories (σ = 0 at the bottom). In each cell `mean` (`error_mean`), `q99` (`error_q99`) and `n.e.` (`non_estimable_area_%`); its color is `error_mean`, sequential colormap, log scale, with the color bar |
| `plot_core.png` | the core sample, as a function of r (1.2, 2, 2.8), k = 4 (15 subsets), σ = 0.05 ms: left `error_mean` and `error_q99`, right the non-estimable area |
| `plot3_error_k_ratio.png` | boxes of the log-Euclidean error E (log) grouped by k = 3…6, in each group one box per scenario of `CORE` (one color per r, legend), on the estimable faces of all the subsets with that k; box q25–q75, line at q50, whiskers q10–q90, area-weighted; in the title Spearman's ρ(E, k) and ρ(E, r) over the same faces, unweighted |
| `plot4_roney_noise_*` | per-face vs Roney, one point per triplet and scenario of `GRID`, 80 points, one Okabe-Ito color and one marker per σ (`σ = 0 ms`, `σ = 0.05 ms`, …) |
| `plot4_roney_ratio_*` | the same on the scenarios of `CORE`, 60 points, one color and marker per r (`r = 1.2`, …) |

Plot 4 is two PNGs per version: `_area`, the non-estimable area of each
method alone, one panel; `_error`, three panels q50, q90 and q99 on the paired
faces, each with its own scale. Roney on the horizontal axis and per-face on
the vertical axis, both linear from 0 with the same limits, with y = x; limits
and centroid over all the points of the version; on the axes only 0, minimum,
centroid (mean of all the points, cross) and maximum, with dotted projections
and without grid.

## Error PLYs

For MeshLab, which colors the faces with a single field, `quality`:

```bash
./estimator/build/ErrorMaps MESH.ply NUMBERS.txt OUT.ply
```

with `MESH.ply` a PLY of the mesh (`error_maps.py` passes
`data/reference/SCENARIO/lat.ply`; only vertices and faces are read) and
`NUMBERS.txt` an `error_maps.txt` written by `error_maps.py`: only numbers, one
row per face in `face_id` order, `estimable error`, with `0 0`
on the non-estimable faces (no `nan`, which `>>` in C++ does not read). The
numbers come from `load()` of `analysis/analysis.py` and are written with
`repr`, which reads back exactly as a double. `ErrorMaps` writes with
`exportAsPly` (`estimator/src/io.cpp`) one text PLY, `OUT.ply` (`format ascii 1.0`,
written by happly, which adds a `comment` line), with the log-Euclidean error
as the per-face `quality`. Header:

```
element vertex 21811
property float x
property float y
property float z
property float nx
property float ny
property float nz
element face N
property list uchar int vertex_indices
property float quality
```

- Vertices: all of them, in the order of `MESH.ply`, `x y z` in mm, with the
  per-vertex normals `nx ny nz`, unit and pointing outward, computed on the
  whole mesh, so a vertex left without faces still has a finite normal.
- Faces: only the estimable ones (`estimable = 1` in `NUMBERS.txt`), in
  `face_id` order; `N` is their number. The triangles
  are those of `MESH.ply`, as they are (outward normals, as in `truth.ply`), so
  MeshLab shows the surface seen from outside. The non-estimable faces are left
  out because MeshLab crashes on a NaN quality ("Colorize by face quality"):
  the holes in the surface are the non-estimable area. A written face with a
  non-finite quality stops `ErrorMaps` with an error naming the file and the
  `face_id`.
- `quality` is the value of `NUMBERS.txt` rounded to float.

These PLYs are only for visualization: they must never be used as input to
any computation (missing faces, float values).

`error_maps.py` chooses the cases, writes for each one `error_maps.txt` next to
its `estimate.ply`, creates the folder of the case
`output/mesh/METHOD_SCENARIO_SUB/` (CASEDIR) and prints a line
`MESH NUMBERS CASEDIR`; `run_all.sh` calls `ErrorMaps MESH NUMBERS CASEDIR/error.ply`.
The cases:

For each scenario of `GRID` (imported from `analysis.py`):

1. per-face, all maps;
2. per-face and Roney on a single triplet, the same at every noise level, so
   that in MeshLab the levels compare on the same maps: the triplet with the
   largest paired area (faces estimable by both methods, area-weighted, areas
   from `load()`) on `painted` without noise.

## Conduction fields for MeshLab

`./estimator/build/FiberView D.ply PREFIX` turns any format D file (ground
truth or estimate) into `PREFIX_ratio.ply` (quality per face = r) and
`PREFIX_fibers.ply` (a thin triangle along each fiber, length proportional to
CV_L); faces with NaN are left out.

`run_all.sh` calls it for the cases of the error maps (`error_maps.py`): for
each one on the estimate, writing `CASEDIR/estimate_{ratio,fibers}.ply`,
and on the ground truth of its scenario, writing
`output/mesh/truth_SCENARIO/truth_{ratio,fibers}.ply` (the four scenarios of
the grid; the same files are written once per case of the scenario).

To compare ground truth and estimate in MeshLab: *Colorize by face quality*
on the two `_ratio.ply`, with Min and Max set by hand and equal for both, and
Percentile Crop 0, so the same colour is the same r.
