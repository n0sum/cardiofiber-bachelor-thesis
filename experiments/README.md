# Experiments

The pipeline, for each scenario (`ALL_SCENARIOS` in `run_all.sh`, built from
`NOISES`, `RATIOS`, `CORE_NOISE` and `CORE_K` at its top; in this order `painted`,
`painted_noise0.05`, `painted_noise0.1`, `painted_noise0.3`,
`painted_r1.2_noise0.05`, `painted_r2_noise0.05`, `painted_r2.8_noise0.05`):

```
simulator/make_truth.py TRUTH            ->  data/reference/TRUTH/truth.ply      ground truth, format D
Simulate truth.ply SITES lat.ply SIGMA   ->  data/reference/SCENARIO/lat.ply     six LAT maps
experiments/subsets.py SCENARIO [K]      ->  output/SCENARIO/confs/sub_….conf     (K = 4 on the core sample)
Estimate output/SCENARIO/confs/*         ->  output/SCENARIO/per_face/sub_…/estimate.{ply,csv}
estimator/roney.py roney.conf SC..       ->  output/SCENARIO/roney/sub_i-j-k/estimate.{ply,csv}   (grid only)
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

The `.conf` files are generated by `subsets.py SCENARIO [K]`, with `write_conf`:
the lines of `template.conf` (the threshold, the same for all) plus `lat_ply`,
`maps`, `output_ply` and `output_csv`. They are not in git: the recipe is the
script. Without K, `subsets.py` writes one `.conf` for each subset of at least
3 maps (the rows of `pacing_sites.csv`; with six maps 20 + 15 + 6 + 1 = 42);
with K, only for each subset of exactly K maps (with six maps and K = 4, 15).
`run_all.sh` calls it without K on the scenarios of the grid and with
`CORE_K` (4, equal to `CORE_K` in `analysis.py`) on those of the core sample.
The `.conf` files go to `output/SCENARIO/confs/`; `subsets.py` also creates the output folder
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
`estimator/README_roney.md`. `run_all.sh` runs it only on the scenarios of the
grid (`GRID_SCENARIOS`: `painted` and `painted_noiseN`); the core sample in
anisotropy is per-face only and has no `roney/` folder.

## The experiment grid

```
output/
  SCENARIO/                           one per scenario
    confs/sub_….conf                  subsets.py, 42 subsets (grid) or 15 with 4 maps (core sample)
    per_face/sub_…/estimate.ply       Estimate, the same subsets
    per_face/sub_…/estimate.csv
    roney/sub_i-j-k/estimate.ply      estimator/roney.py, 20 triplets (scenarios of the grid only)
    roney/sub_i-j-k/estimate.csv
  figures/                            analysis/analysis.py, tables and PNGs
  mesh/                               all the PLYs for MeshLab (ErrorMaps, FiberView), three folders, below
```

In the folders of the error-map cases (`error_maps.py`) there is also
`error_maps.txt`. `output/mesh/`, with TRIPLET and the subsets SUB chosen by
`error_maps.py` (below):

```
output/mesh/
  roney_comparison/   per_face_S_TRIPLET_error.ply  roney_S_TRIPLET_error.ply     for each scenario S of the grid (8)
  worst_subsets/      per_face_painted_noise0.05_kK_SUB_error.ply                 K = 3, 4, 5, 6 (4)
  ellipses/           truth_painted_ellipses.ply
                      per_face_painted_noise0.05_k3_SUB_ellipses.ply              SUB of worst_subsets k3
                      per_face_painted_noise0.05_k6_sub_0-1-2-3-4-5_ellipses.ply
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
`run_all.sh`, for `per_face.csv`. `comparison.csv` and plot 4 are on the
scenarios of `GRID` only, the only ones with Roney.

`load(scenario, method, sub)` reads `data/reference/SCENARIO/truth.ply` (once per
scenario), `estimate.ply` and `estimate.csv` of the run, and returns per face:

| column      | meaning |
|-------------|---------|
| `area`      | area of the face from the vertices of the ground truth, `|(b − a) × (c − a)| / 2`, mm² |
| `estimable` | μ of the estimate finite |
| `cond` / `k_used` | from `estimate.csv` |
| `error`     | log-Euclidean distance between estimated and true tensor, below |
| `angle`, `mean_speed`, `ratio` | the three errors of the core sample, below |

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

For the core sample, where r changes and E weights the fiber with ln r,
`load()` also gives three errors that separate the fiber from the speeds, with
v = √(CV_L·CV_T) = e^μ the mean speed and r = CV_L/CV_T = e^(2α):

| column       | meaning |
|--------------|---------|
| `angle`      | angle between estimated and true fiber, degrees of arccos(min(1, \|f̂·f\|)), unoriented |
| `mean_speed` | 100 · (v̂ − v) / v, % |
| `ratio`      | 100 · (r̂ − r) / r, % |

The errors are NaN on the non-estimable faces. Statistics are area-weighted:
q50, q90 and q99 (`error_q50/q90/q99` in the tables) are the values below
which 50, 90 and 99% of the area lies.

- `per_face.csv`: one row per scenario and subset of the per-face
  estimator (on the core sample only the 15 subsets with 4 maps): non-estimable area, errors on the estimable faces of the
  subset, `cond_q50/q90/q99` on the faces with computed cond.
- `comparison.csv`: one row per scenario of `GRID`, triplet and method (`per_face`,
  `roney`), with the non-estimable area of the method alone, the paired area
  (faces estimable by both methods with that triplet) and the errors on the
  paired faces. For the Roney rows the non-estimable area is also split in
  two: `fit_area_%` (`k_used < 3`, lost to the planar fit) and `conic_area_%`
  (`k_used = 3` and estimate NaN: the conic is not an ellipse or the system is
  nearly singular); both empty on the per-face rows. On the paired faces
  `error_mean` = Σ area·E / Σ area and `error_q50/q90/q99`. The rows with `sub` =
  `pooled` gather the paired faces of the 20 triplets of the scenario; their
  areas are means over the triplets.
- `grid.csv`: one row per (σ, k) of the scenarios of `GRID`, columns
  `scenario, noise_ms, k, subsets, non_estimable_area_%, error_mean,
  error_q50, error_q90, error_q99`.
- `core.csv`: one row per scenario of `CORE`, at k = `CORE_K`, columns
  `scenario, r, k, subsets, non_estimable_area_%` and, for each of `angle`,
  `mean_speed` and `ratio`, `_mean`, `_q50`, `_q90`, `_q99` of the absolute
  value (e.g. `angle_mean`, `mean_speed_q99`, `ratio_q50`). No E.

A cell (scenario, k) of `grid.csv` and `core.csv` (`cell()`) takes all the
subsets with k maps of the scenario (20 triplets, 15 quadruplets, 6
quintuplets, 1 sextuplet) and concatenates their DataFrames from `load()`:
each face appears once per subset. For each error X (E in `grid.csv`; angle,
mean speed and ratio in `core.csv`) `X_mean` = Σ area·|X| / Σ area and
`X_q50/q90/q99` (area-weighted quantiles of |X|) are on the estimable rows;
`non_estimable_area_%` = area of the non-estimable rows / area of all the
rows × 100, i.e. the mean over the subsets of the non-estimable area of each;
`subsets` = number of subsets of the cell. The cond (`cond_q50/q90/q99` of
`per_face.csv`, from `estimate.csv`) is in no figure.

| file | content |
|------|---------|
| `plot_grid.png` | the grid: a 4 × 4 colored table, k = 3…6 on the horizontal axis, σ = 0, 0.05, 0.1, 0.3 ms on the vertical axis as equispaced categories (σ = 0 at the bottom). In each cell `mean` (`error_mean`), `q99` (`error_q99`) and `n.e.` (`non_estimable_area_%`); its color is `error_mean`, sequential colormap, log scale, with the color bar |
| `plot_core.png` | the core sample, as a function of r (1.2, 2, 2.8), k = 4 (15 subsets), σ = 0.05 ms, 2 × 2 panels: fiber angle [°], relative error on the mean speed \|v̂ − v\|/v [%], relative error on the ratio \|r̂ − r\|/r [%], each with mean (solid) and q99 (dashed), and the non-estimable area |
| `plot_triplets_area.png` | per-face, the 20 triplets of each scenario of `GRID`: estimable area (100 − `non_estimable_area_%` of `per_face.csv`, k = 3) as a function of σ, as equispaced categories (0, 0.05, 0.1, 0.3 ms), linear up to 100%; one point per triplet and level, the points of a triplet joined by a thin gray line, staggered horizontally by an offset that depends only on the triplet |
| `plot4_roney_noise_*` | per-face vs Roney, one point per triplet and scenario of `GRID`, 80 points, one Okabe-Ito color and one marker per σ (`σ = 0 ms`, `σ = 0.05 ms`, …) |

Plot 4 is two PNGs: `_area`, the non-estimable area of each
method alone, one panel; `_error`, two panels on the paired faces (the faces
estimated by both methods), mean (`error_mean`) and q99 (`error_q99`), each
with its own scale. Roney on the horizontal axis and per-face on the vertical
axis, both linear from 0 with the same limits, with y = x; limits over all the
points; on the axes only the minimum and the maximum of the cloud (no 0 tick),
with dotted projections and without grid. Legend: the four σ and y = x.

## PLYs for MeshLab

Every PLY in `output/mesh/` is text (`format ascii 1.0`, written by happly,
which adds a `comment` line) and carries its own color, one RGB per face,
written by `writeColoredPly` (`estimator/src/io.cpp`), shared by `ErrorMaps`
and `FiberView`. Header:

```
element vertex NV
property float x
property float y
property float z
element face NF
property list uchar int vertex_indices
property uchar red
property uchar green
property uchar blue
```

Vertices `x y z` in mm, no normals; the faces of the mesh are those of the
input, in `face_id` order and as they are (outward normals), so MeshLab shows
the surface seen from outside. In MeshLab: *Render → Color → Per Face*.
The colors (constants at the top of `ErrorMaps.cpp` and `FiberView.cpp`),
rounded to integers 0-255:

- gray (160, 160, 160): non-estimable face (error maps) or mesh face (ellipses);
- error E: t = min(E / E_MAX, 1), red = 255 min(1, 2t), green = 255 min(1, 2(1 − t)),
  blue = 0; green = 0, yellow = E_MAX / 2 ≈ 0.98, red = E_MAX ≈ 1.96 and above, with
  E_MAX = 2√2 ln 2, the error of a fiber rotated by 90° with r = 2;
- mean speed v = e^μ: s = clamp((v − 0.5) / (1.0 − 0.5), 0, 1), red = 255 s,
  green = 0, blue = 255 (1 − s); blue = 0.5 mm/ms and below, red = 1.0 mm/ms and
  above, never through gray.

These PLYs are only for visualization: they must never be used as input to
any computation (colors, float coordinates).

### Error maps

```bash
./estimator/build/ErrorMaps MESH.ply NUMBERS.txt OUT.ply
```

with `MESH.ply` a PLY of the mesh (`error_maps.py` passes
`data/reference/SCENARIO/lat.ply`; only vertices and faces are read) and
`NUMBERS.txt` an `error_maps.txt` written by `error_maps.py`: only numbers, one
row per face in `face_id` order, `estimable error`, with `0 0`
on the non-estimable faces (no `nan`, which `>>` in C++ does not read). The
numbers come from `load()` of `analysis/analysis.py` and are written with
`repr`, which reads back exactly as a double. `OUT.ply` has all the vertices
(21811) and all the faces (43270) of `MESH.ply`, no holes: gray where
`estimable = 0`, otherwise the color of the error. An estimable face with a
non-finite error stops `ErrorMaps` with an error naming the file and the
`face_id`.

`error_maps.py` (no arguments) chooses the cases, writes for each error map
`error_maps.txt` next to its `estimate.ply`, creates the three folders of
`output/mesh/` and prints one line per case, the kind first:
`error MESH NUMBERS OUT` or `ellipses D OUT`, with OUT the final path of the
PLY; `run_all.sh` calls `ErrorMaps MESH NUMBERS OUT` or `FiberView D OUT`.
The cases (`GRID` imported from `analysis.py`):

1. `roney_comparison/`: for each scenario of `GRID`, per-face and Roney on a
   single triplet, the same at every noise level, so that the levels compare
   on the same maps: the triplet with the largest paired area (faces estimable
   by both methods, area-weighted, areas from `load()`) on `painted` without
   noise;
2. `worst_subsets/`: on the scenario of `GRID` with σ = `WORST_NOISE` = 0.05 ms,
   per-face, for each k = 3…6 the subset with k maps with the highest mean
   error, sum(area · error) / sum(area) on its estimable faces (with k = 6 the
   only one, `sub_0-1-2-3-4-5`);
3. `ellipses/`: the ground truth `data/reference/painted/truth.ply`, and on
   `painted_noise0.05` the per-face estimates of the triplet chosen for k = 3
   in `worst_subsets` and of `sub_0-1-2-3-4-5`.

### Conduction fields as ellipses

```bash
./estimator/build/FiberView D.ply OUT.ply
```

turns a format D file (ground truth or estimate) into one PLY with, in the
same file, the mesh of `D.ply` with all its faces gray and, for each face with
a value (`fx fy fz mu alpha` all finite), an ellipse: a fan of
`N_SEGMENTS` = 16 triangles around the centroid c of the face, lifted by
`H` = 0.05 mm along the outward normal n. With f the fiber projected on the
plane of the face and normalized and t = n × f, the vertices are c + H n and,
for i = 0…15, θ = 2π i / 16, c + `RHO` cos θ f + (`RHO` / r) sin θ t + H n,
with `RHO` = 0.23 mm and r = e^(2α); the triangles (center, i, i + 1 mod 16)
have normal +n. The long axis is fixed: every ellipse lies inside the circle
of radius `RHO`, a circle with α = 0, a needle 2 `RHO` long with large α; the
shape shows fiber and anisotropy, the color the mean speed v = e^μ (all 16
triangles). Faces without a value: only gray, no ellipse. Vertices
21811 + 17 · (faces with a value), faces 43270 + 16 · (faces with a value).
