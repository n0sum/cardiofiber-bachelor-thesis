# Differential estimation of atrial conduction anisotropy

Bachelor's thesis in geometry processing. On a left atrium with known fibers
and conduction velocities (CV), six activation maps (LAT) are simulated with an
eikonal solver; from the maps alone the conduction tensor (fiber direction,
CV_L, CV_T) is estimated per face with the per-face gradient, and the estimate
is validated against the known ground truth, also in comparison with method 2
of Roney et al. 2019 (*Computers in Biology and Medicine* 104, 278-290).

```
run_all.sh                     rebuilds data/reference/ and output/ from data/
environment.yml                conda environment "cardiofiber"
data/
  external/biatrial_model_completeFit.h5      SSM (included)
  augmenta/LA_endo_with_fiber.{pts,elem,lon}  LA with fibers, AugmentA output (included)
simulator/
  make_truth.py                ground truth of a scenario in format D
  painted_{alpha,mu}.ply       alpha and mu painted in MeshLab, scenario painted (included)
  stage1/                      LA from the SSM (optional, see below), pacing sites
cardiofiber/                   estimation package (C++ library, glm only): mesh, per-vertex LAT maps
                               -> per-face gradient -> per-face tensor; cardiofiber/README.md
estimator/                     project around the package: files and executables
  src/                         Simulate (C++): eikonal LAT maps from a format D field;
                               Estimate (C++): experiments from the .conf files, estimate in format D;
                               ErrorMaps (C++): error PLYs for MeshLab;
                               FiberView (C++): a format D field as PLYs for MeshLab
  include/happly.h             vendored PLY reader/writer (THIRD_PARTY.md)
  roney.py                     method 2 of Roney et al. 2019, estimate in format D
analysis/analysis.py           errors against the ground truth, area-weighted statistics,
                               figures, per_face.csv and comparison.csv
experiments/                   template.conf and roney.conf, subsets.py, error_maps.py;
                               formats in experiments/README.md
```

## Requirements

- Linux; developed on Ubuntu 26.04.1.
- conda, with the environment of `environment.yml` (Python 3.13.15, numpy, scipy,
  pandas, matplotlib, h5py, meshio, pyvista 0.48.4, vtk 9.6.2: exact
  versions in the file).
- CMake >= 3.20 (tested 4.2.3), a C++20 compiler (tested GCC 15.2.0) and
  glm, header-only (on Ubuntu `libglm-dev`, tested 1.0.1).
- Docker, only to optionally regenerate the AugmentA output.

## Included data

These are the only inputs that are not regenerated; everything else in `data/` and
`output/` is written by `run_all.sh`.

**SSM** `data/external/biatrial_model_completeFit.h5` (23.9 MB). Bi-atrial
statistical shape model from KIT-IBT, Zenodo
<https://zenodo.org/records/4309958>, CC-BY 4.0 license.
md5 `96d8bf64bc943dcba71e3702a552fa75`,
sha256 `d4b3a90f2b9c6db79bb225954a4e5eeaed59f44916e5a93935dde20458da63b8`.
Citation: C. Nagel, S. Schuler, O. Dössel, A. Loewe, *A bi-atrial
statistical shape model for large-scale in silico studies of human atria:
model development and application to ECG simulations*, Medical Image
Analysis 74 (2021) 102210.

**AugmentA output** `data/augmenta/LA_endo_with_fiber.{pts,elem,lon}`
(3.7 MB): the mean left atrium of the SSM, 21811 vertices and 43270
triangles, with fibers and per-triangle region tags. Generated from the SSM with
`make_mean.py`, `extract_la.py` and then AugmentA with the flags
of the section *Regenerating the LA and the fibers*. 
AugmentA license: Academic Public License. Citation:
L. Azzolin et al., *AugmentA: patient-specific augmented atrial model
generation tool*, Computerized Medical Imaging and Graphics 108 (2023)
102265.

**Painted maps** `simulator/painted_alpha.ply` and `simulator/painted_mu.ply`:
α and μ of the `painted` scenario, painted by hand in MeshLab in grayscale
on the AugmentA mesh with flipped faces (21811 vertices, per-vertex colors).
`make_truth.py` reads only their colors, matched by vertex index; positions,
faces and fibers come from `data/augmenta/`.

## Reproducing the results

```bash
git clone https://github.com/n0sum/cardiofiber
cd cardiofiber
conda env create -f environment.yml
conda activate cardiofiber
```

Then from the root, with the `cardiofiber` environment active:

```bash
./run_all.sh
```

The script builds `estimator/` in `estimator/build/`, then for each scenario
(the list `SCENARIOS` at its top) writes the ground truth with
`make_truth.py` and simulates the six LAT maps with `Simulate` (one pacing
site per row of `simulator/stage1/pacing_sites.csv`). It then deletes
`output/` and computes a single experiment grid: for each scenario the
per-face estimator on all 42 subsets of at least 3 maps and Roney's method on
the 20 triplets. The two estimators read only the LAT maps; the errors are
computed afterwards by `analysis/analysis.py`, against the ground truth.

Scenarios (`simulator/make_truth.py`):

| scenario  | ground truth |
|-----------|--------------|
| `uniform` | AugmentA fibers, CV_L = 1.05, CV_T = 0.525 mm/ms on every face (r = 2) |
| `painted` | AugmentA fibers; μ and α vary continuously, painted by hand (`simulator/painted_{alpha,mu}.ply`) |
| `r1.6`, `r2.4` | sweep in r, with `uniform` as r = 2: as `uniform` (same fibers, same μ), α = ½ ln r on every face with r = 1.6, 2.4 |

In `painted`, per vertex g = mean of R, G, B / 255 (0 black, 1 white; the
alpha channel is ignored), then α = ½ ln 1.6 + g (½ ln 2.4 − ½ ln 1.6)
(r from 1.6, black, to 2.4, white) and μ = −0.6 + 0.6 g (e^μ from 0.55,
black, to 1.0 mm/ms, white); per face the mean of the values of its three
vertices.

Outputs:

| folder                     | content |
|----------------------------|---------|
| `data/reference/SCENARIO/` | `truth.ply` (ground truth, format D), `lat.ply` (six LAT maps) |
| `output/SCENARIO/`         | `confs/sub_*.conf`; `per_face/sub_*/estimate.{ply,csv}` (Estimate); `roney/sub_i-j-k/estimate.{ply,csv}` (`roney.py`); `error_maps.txt` (`error_maps.py`) next to the `estimate.ply` of the error-PLY cases |
| `output/figures/`          | `per_face.csv`, `comparison.csv`, `cond_k_ratio.csv` and the plots (`analysis.py`) |
| `output/mesh/`             | all the PLYs for MeshLab, one folder per case and one per ground truth (`ErrorMaps`, `FiberView`) |

`output/mesh/`, with the cases chosen today by `error_maps.py` (the triplet is the one
of `painted` with the largest paired area):

```
output/mesh/
  per_face_uniform_sub_0-1-2-3-4-5/   error_angle.ply  error_mean_cv.ply  error_ratio.ply  estimate_ratio.ply  estimate_fibers.ply
  per_face_painted_sub_0-1-2-3-4-5/   (the same five)
  per_face_painted_sub_1-4-5/         (the same five)
  roney_painted_sub_1-4-5/            (the same five)
  truth_uniform/                      truth_ratio.ply  truth_fibers.ply
  truth_painted/                      truth_ratio.ply  truth_fibers.ply
```

To compare a `truth_ratio.ply` with an `estimate_ratio.ply` in MeshLab:
*Colorize by face quality* with Min and Max set by hand, equal for both, and
Percentile Crop 0.

The numbers are in `output/figures/`: `per_face.csv` for the per-face
estimator (one row per scenario and subset), `comparison.csv` for the
comparison with Roney on the map triplets. Errors are area-weighted over the
estimable faces: the angle, in degrees, between the estimated and the true
fiber; for the mean speed e^μ = √(CV_L·CV_T) and for the ratio r = CV_L / CV_T
the difference from the true value, as a percentage of the true value. Figures and columns are
described in `experiments/README.md`.

## Units and formats

- Lengths in mm, times in ms, velocities in mm/ms, everywhere: SSM, AugmentA
  output, ground truth, LAT maps and estimates.
- `.pts`: first line the number of vertices, then `x y z` per line.
- `.elem`: first line the number of triangles, then `Tr i j k tag`, 0-based indices,
  the tag is the AugmentA region (1-7, not used). AugmentA triangles have inward
  normals; `make_truth.py` flips them once, and every file after it has
  outward normals.
- `.lon`: first line `2`, then per triangle fiber and sheet, three components each.
- Format D (`truth.ply`, `estimate.ply`): binary PLY, per vertex `x y z`, per face
  the triangle, the unit fiber `fx fy fz` and `mu`, `alpha` with
  CV_L = e^(mu+alpha), CV_T = e^(mu-alpha); NaN on faces without a value.
- `lat.ply`: binary PLY, per vertex `x y z lat_0 … lat_5` (mm, ms), faces as in
  `truth.ply`; `lat_k` from site k.
- `estimate.csv`, the `.conf` files, `error_maps.txt` and the PLYs for MeshLab:
  `experiments/README.md`.

## Regenerating the LA and the fibers from the SSM (optional)

`run_all.sh` does not require it.

```bash
python simulator/stage1/make_mean.py      # -> data/stage1/ssm_mean_surface.vtk
python simulator/stage1/extract_la.py     # -> data/stage1/ssm_mean_LA.vtk, 21811 vertices / 43270 triangles
cp simulator/stage1/ssm_mean_LA_apexes.csv data/stage1/
```

`extract_la.py` rewrites `simulator/stage1/la_vertex_mask.npy`: `git status`
must stay clean. Then AugmentA (<https://github.com/KIT-IBT/AugmentA>), with
the `docker/Docker-compose.yml` of its repository and the data volume pointed
at `data/stage1/`:

```yaml
    volumes:
      - /absolute/path/to/cardiofiber/data/stage1:/headless/data
      - ../:/headless/augmenta
```

```bash
docker build -f Dockerfile-x64 -t augmenta_container_image .   # from the docker/ folder of AugmentA
docker compose -f Docker-compose.yml up
docker exec -u root CONTAINER pip install ruamel.yaml           # CONTAINER from docker ps
```

In the container desktop (VNC `localhost:5901`, or `localhost:6901/?password`):

```bash
cd /headless/augmenta
python main.py --mesh /headless/data/ssm_mean_LA.vtk --closed_surface 0 --open_orifices 0 --atrium LA --SSM_fitting 0 --resample_input 0 --normals_outside 1 --add_bridges 0 --MRI 0 --ofmt vtk --apex-file /headless/data/ssm_mean_LA_apexes.csv
```

`--resample_input 0` is essential: it keeps the mesh in the SSM indices.
`--apex-file` replaces the interactive selection of the apex
of the appendage done in the original run (2026-08-20): the point closest
to the coordinates in the file is vertex 1405, the same one. The actual config
of the original run is `simulator/stage1/augmenta_config.yaml`. The fibers
come out in `data/stage1/ssm_mean_LA_fibers/result_LA/LA_endo_with_fiber.{pts,elem,lon}`.

## License

The code of the repository is under AGPL-3.0-or-later (see [`LICENSE`](LICENSE)),
including the estimation package `cardiofiber/`.
The following are exceptions, and remain under their own licenses:

- `estimator/include/happly.h`: MIT, license in the file header (see
  [`estimator/THIRD_PARTY.md`](estimator/THIRD_PARTY.md));
- the SSM in `data/external/`: CC-BY 4.0, Nagel et al. 2021;
- the AugmentA output in `data/augmenta/`: generated with AugmentA, Academic
  Public License.
