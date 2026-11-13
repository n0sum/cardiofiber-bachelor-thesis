#include <fstream>
#include <string>
#include <vector>
#include <cardiofiber/Mesh.h>
#include "io.h"

// Usage, from the repository root:
//   ./estimator/build/ErrorMaps MESH.ply NUMBERS.txt PREFIX
// MESH.ply: a PLY of the mesh (vertices and faces; e.g. a lat.ply), faces used as they are, nothing is flipped.
// NUMBERS.txt: one row per face, "estimable angle mean_cv ratio" (experiments/error_maps.py).
// Writes PREFIX_{angle,mean_cv,ratio}.ply. The error is the quality of the face; only
// the faces with estimable = 1 are written, the same in the three files: a NaN quality crashes MeshLab.
// The vertices are all written, so the face indices do not change.

int main(int argc, char** argv) {
    const MultiLatPly in = readMultiLatPly(argv[1], {});
    const cardiofiber::Mesh mesh(in.vertices, in.faces, {});

    std::vector<bool> estimable;
    std::vector<double> angle, meanCv, ratio;
    std::ifstream numbers(argv[2]);
    int e;
    double a, m, r;
    while (numbers >> e >> a >> m >> r) {
        estimable.push_back(e == 1);
        angle.push_back(a);
        meanCv.push_back(m);
        ratio.push_back(r);
    }

    const std::string prefix = argv[3];
    exportAsPly(prefix + "_angle.ply", mesh, angle, estimable);
    exportAsPly(prefix + "_mean_cv.ply", mesh, meanCv, estimable);
    exportAsPly(prefix + "_ratio.ply", mesh, ratio, estimable);
    return 0;
}
