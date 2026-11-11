# cardiofiber

Face-by-face estimation of the conduction tensor on a triangle mesh
from several activation time (LAT) maps, one LAT per vertex per
map:

    mesh + per-vertex LAT  ->  per-face gradient  ->  per-face tensor

For each face the output is the fiber direction (`fiber_dir`, unit, tangent
to the face, unoriented), the velocity along the fiber `CV_L` and the
transverse one `CV_T`, the condition number `cond` of the system and the flag
`estimable`.

It is a C++20 static library that depends only on [glm](https://github.com/g-truc/glm):
it takes in-memory lists and neither reads nor writes files. Everything is in the namespace
`cardiofiber`.

## Usage

```cpp
#include <cardiofiber/Mesh.h>

// vertices: std::vector<glm::dvec3>; faces: std::vector<cardiofiber::Face> (indices a, b, c);
// lats: std::vector<std::vector<double>>, lats[map][vertex]
const cardiofiber::Mesh m(vertices, faces, lats);
const std::vector<std::vector<glm::dvec3>> gradients = m.computeGradients();      // gradients[map][face]
// threshold on sigma_min; one FaceAnisotropy per face
const std::vector<cardiofiber::FaceAnisotropy> est = m.estimateAnisotropy(gradients, 1e-3);
```

With CMake:

```cmake
add_subdirectory(cardiofiber)
target_link_libraries(my_program PRIVATE cardiofiber)
```

The constructor copies the lists and throws `std::invalid_argument` if a map
does not have exactly one LAT per vertex; `estimateAnisotropy` does the same if a
map does not have exactly one gradient per face.

## Units

Velocities come out in coordinate units / LAT units: with
coordinates in mm and LAT in ms, as in the rest of the repository, in mm/ms.

## Gradients given from outside

`estimateAnisotropy` accepts any gradients, not only those from
`computeGradients()`. They must be given in 3D and tangent to the face: the estimate
projects them onto the local basis of the face and the component along the normal is
discarded without warning.

## Non-estimable faces

A face has `estimable = false`, and the other fields not meaningful, if:

- it has fewer than 3 maps (underdetermined system);
- sigma_min, the smallest singular value of the system, does not exceed the
  threshold `minSingularValue`;
- lambda2 = CV_T^2 <= 1e-8 lambda1, with lambda1 = CV_L^2: a test of numerical zero
  (CV_T / CV_L > 1e-4).

The threshold on sigma_min is absolute, not relative: the rows of the system have
the units of the gradient squared, so the right threshold depends on the
scale of the velocities (and on the chosen units).

## Limitations

Validated only on exact and dense LATs (one per vertex, from simulation), without
noise.

## License

AGPL-3.0-or-later, like the rest of the repository (see [`../LICENSE`](../LICENSE)).
