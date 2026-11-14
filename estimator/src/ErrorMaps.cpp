#include <fstream>
#include <vector>
#include <cardiofiber/Mesh.h>
#include "io.h"

// Usage, from the repository root:
//   ./estimator/build/ErrorMaps MESH.ply NUMBERS.txt OUT.ply
// MESH.ply: a PLY of the mesh (vertices and faces; e.g. a lat.ply), faces used as they are, nothing is flipped.
// NUMBERS.txt: one row per face, "estimable error" (experiments/error_maps.py), the log-Euclidean
// distance between estimated and true tensor.
// Writes OUT.ply. The error is the quality of the face; only the faces with estimable = 1 are
// written: a NaN quality crashes MeshLab. The vertices are all written, so the face indices do not change.

int main(int argc, char** argv) {
    const MultiLatPly in = readMultiLatPly(argv[1], {});
    const cardiofiber::Mesh mesh(in.vertices, in.faces, {});

    std::vector<bool> estimable;
    std::vector<double> error;
    std::ifstream numbers(argv[2]);
    int e;
    double v;
    while (numbers >> e >> v) {
        estimable.push_back(e == 1);
        error.push_back(v);
    }

    exportAsPly(argv[3], mesh, error, estimable);
    return 0;
}
