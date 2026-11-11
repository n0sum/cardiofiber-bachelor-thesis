# Differential estimation of atrial conduction anisotropy

Bachelor's thesis in geometry processing. On a left atrium with known fibers
and conduction velocities (CV), six activation maps (LAT) are simulated;
from these the conduction tensor (fiber direction, CV_L, CV_T) is estimated
per face with the per-face gradient, and the estimate is validated against
the known ground truth, also in comparison with method 2 of Roney et al. 2019
(*Computers in Biology and Medicine* 104, 278-290).

```
run_all.sh                     rebuilds data/reference/ and output/ from data/
environment.yml                conda environment "cardiofiber"
data/
  external/biatrial_model_completeFit.h5      SSM (included)
  augmenta/LA_endo_with_fiber.{pts,elem,lon}  LA with fibers, AugmentA output (included)
simulator/
  stage1/                      LA from the SSM (optional, see below), pacing sites
    opencarp/                  mesh in um, per-zone tags, eikonal simulations
  convert/                     from openCARP to PLY: LAT maps and per-face ground truth
cardiofiber/                   estimation package (C++ library, glm only): mesh, per-vertex LAT maps
                               -> per-face gradient -> per-face tensor; cardiofiber/README.md
estimator/                     project around the package: files, ground truth, errors
  src/                         Estimate (C++): experiments from the .conf files, per-face errors.csv;
                               ErrorMaps (C++): error PLYs for MeshLab
  include/happly.h             vendored PLY reader (THIRD_PARTY.md)
  roney.py                     method 2 of Roney et al. 2019
analysis/analysis.py           area-weighted error statistics, figures, per_face.csv and comparison.csv
experiments/                   template.conf and roney.conf, per-zone CV (region_cv/),
                               subsets.py, error_maps.py; formats in experiments/README.md
```

## Requirements

- Linux; developed on Ubuntu 26.04.1.
- conda, with the environment of `environment.yml` (Python 3.13.15, numpy, scipy,
  pandas, matplotlib, h5py, meshio, pyvista 0.48.4, vtk 9.6.2: exact
  versions in the file).
- openCARP v19.0, from the standard `.deb` package (not ubuntu20):
  <https://opencarp.org/download/installation>.
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

## Reproducing the results

```bash
git clone https://github.com/n0sum/cardiofiber
cd cardiofiber
conda env create -f environment.yml
conda activate cardiofiber
```

Install openCARP v19.0 (see *Requirements*), then from the root, with
the `cardiofiber` environment active:

```bash
./run_all.sh
```

The script converts the mesh to um, writes the scenario tags, runs the
eikonal simulations (six pacing sites, `simulator/stage1/pacing_sites.csv`),
builds `estimator/` in `estimator/build/`, deletes `output/` and computes
a single experiment grid: for each scenario the per-face estimator on
all 42 subsets of at least 3 maps and Roney's method on the 20
triplets. Outputs in `output/`:

| folder                                     | content                                            |
|--------------------------------------------|----------------------------------------------------|
| `aniso_zone_N/`                            | `confs/sub_*.conf`; `per_face/sub_*/errors.csv` (Estimate); `roney/sub_i-j-k/errors.csv` (`roney.py`); `error_maps.txt` next to the `errors.csv` of the error-PLY cases (`error_maps.py`) |
| `figures/`                                 | `per_face.csv`, `comparison.csv` and 19 plots (`analysis.py`) |
| `mesh/`                                    | 12 error PLYs, one per error and case, for MeshLab (`error_maps.py`) |

Scenarios: `aniso_zone_N`, the mesh divided into N random zones (seed 0), each zone
with one of the 25 combinations of CV_L and r = CV_L/CV_T of
`region_cv/aniso_random.csv`, N = 1, 8, 70, 600: mean zone diameter
d = sqrt(4 A / (pi N)) ~ 123 / sqrt(N) mm, i.e. 123, 44, 15 and 5.0 mm, equal
steps on a log scale. With N = 1 the CVs are uniform (tag 16 came out:
CV_L 1.155, CV_T 0.722 mm/ms, r = 1.6); with N = 600 the zones are as large as
the disc of Roney's method (diameter 5 mm), which assumes uniform CVs
within it: with smaller zones the comparison would only measure that the
assumption is violated. Details in `experiments/README.md`.

## Expected results

Errors are area-weighted over the estimable faces. For the fiber direction, the error is the angle, in degrees, between the estimated and the true fiber. For CV_L and for the ratio r, the error is the difference from the true value, taken without sign and expressed as a percentage of the true value.

The tables are in `output/figures/`: `per_face.csv` for the per-face estimator
(one row per scenario and subset), `comparison.csv` for the comparison with
Roney.

**aniso_zone_1**, six maps (`per_face.csv`, row `sub_0-1-2-3-4-5`):

| non-estimable area | angle q50 | angle q90 | CV_L q50 | CV_L q90 | r q50 | r q90 |
|--------|--------|-------|--------|--------|--------|---------|
| 0.13 % | 0.907° | 4.27° | 0.67 % | 3.32 % | 1.30 % | 6.90 %  |

**Map subsets** on `aniso_zone_1` (`per_face.csv`, errors on the
estimable faces of each subset): with 3 maps, over the 20 triplets, non-estimable
area 4.0-9.7 % (median 5.3 %), angle q50 1.17-1.76° (median
1.37°), angle q90 9.2-15.0° (median 10.6°); with 6 maps the values of
`aniso_zone_1`.

**Zones**, six maps (`per_face.csv`, rows `sub_0-1-2-3-4-5`), d = mean zone diameter:

| N   | d [mm] | non-estimable area | angle q50  | angle q90  | CV_L q50 | r q50  |
|-----|--------|--------------------|------------|------------|----------|--------|
| 600 | 5.03   | 1.29 %             | 1.78°      | 9.50°      | 1.99 %   | 4.69 % |
| 70  | 14.7   | 0.57 %             | 1.22°      | 6.50°      | 1.34 %   | 2.99 % |
| 8   | 43.6   | 0.42 %             | 0.882°     | 4.70°      | 1.00 %   | 2.16 % |
| 1   | 123    | 0.13 %             | 0.907°     | 4.27°      | 0.67 %   | 1.30 % |

**Roney vs per-face estimator**, map triplets

| N   | paired area   | method   | non-estimable area | angle q50  | angle q90  | CV_L q50 | r q50  |
|-----|---------------|----------|--------------------|------------|------------|----------|--------|
| 1   | 75.7 %        | Roney    | 22.5 %             | 6.65°      | 21.2°      | 12.9 %   | 10.9 % |
|     |               | per-face | 5.8 %              | 1.31°      | 9.08°      | 0.82 %   | 2.13 % |
| 8   | 55.0 %        | Roney    | 42.4 %             | 9.16°      | 29.1°      | 13.2 %   | 14.6 % |
|     |               | per-face | 9.8 %              | 1.34°      | 9.99°      | 1.08 %   | 3.24 % |
| 70  | 54.8 %        | Roney    | 41.8 %             | 8.75°      | 30.2°      | 13.9 %   | 15.4 % |
|     |               | per-face | 12.6 %             | 1.63°      | 12.4°      | 1.27 %   | 3.84 % |
| 600 | 52.3 %        | Roney    | 42.6 %             | 8.69°      | 29.5°      | 15.0 %   | 17.8 % |
|     |               | per-face | 16.3 %             | 2.38°      | 16.9°      | 1.89 %   | 6.40 % |

## Units and formats

- Mesh and SSM in mm; openCARP wants um. Only lengths change: mm_to_um.py multiplies the coordinates by 1000 and run_sims.py multiplies the velocities by 1000 (mm/ms → um/ms). Time stays in ms, so the LAT that openCARP returns is used as is, and the estimated CVs come out in mm/ms. Unit errors are silent.
- `.pts`: first line the number of vertices, then `x y z` per line.
- `.elem`: first line the number of triangles, then `Tr i j k tag`, 0-based indices. In
  `data/augmenta/` the tag is the AugmentA region (1-7, not used); in the scenario meshes
  `data/reference/aniso_zone_N/` it is the velocity zone, a row of `aniso_random.csv`.
  AugmentA flips the orientation of the triangles with respect to the SSM; the reference
  connectivity is always `LA_endo_with_fiber.elem`.
- `.lon`: first line `2`, then per triangle fiber and sheet, three components each.
- `experiments/region_cv/aniso_random.csv`: `region_tag,CV_L,CV_T` in mm/ms;
  `run_sims.py` converts to mm/s for openCARP.
- `data/reference/LA_LAT_NAME.ply`: binary PLY, per vertex `x y z` (mm) and
  `lat_0 ... lat_5` (ms), `lat_k` from site k. `LA_GT_NAME.ply`: per face `fx fy fz` (fiber),
  `CV_L CV_T` (mm/ms), `elem_tag`.
- `output/aniso_zone_N/{per_face,roney}/sub_*/errors.csv`: one row per face; columns in `experiments/README.md`.

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
