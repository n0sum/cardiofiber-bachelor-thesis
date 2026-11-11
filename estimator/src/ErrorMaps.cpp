#include <fstream>
#include <string>
#include <vector>
#include <cardiofiber/Mesh.h>
#include "io.h"

// Usage, from the repository root:
// The error is the quality of the face; only
// the faces with estimable = 1 are written, the same in the three files: a NaN quality crashes MeshLab.
// The vertices are all written, so the face indices do not change.

int main(int argc, char** argv) {
    const MultiLatPly in = readMultiLatPly(argv[1], {});
    std::vector<cardiofiber::Face> reversedFaces(in.faces.size());
    for (size_t i = 0; i < in.faces.size(); ++i)
        reversedFaces[i] = {in.faces[i].a, in.faces[i].c, in.faces[i].b};
    const cardiofiber::Mesh mesh(in.vertices, reversedFaces, {});

    std::vector<bool> estimable;
    std::vector<double> angle, cvL, ratio;
    std::ifstream numbers(argv[2]);
    int e;
    double a, l, r;
    while (numbers >> e >> a >> l >> r) {
        estimable.push_back(e == 1);
        angle.push_back(a);
        cvL.push_back(l);
        ratio.push_back(r);
    }

    const std::string prefix = argv[3];
    exportAsPly(prefix + "_angle.ply", mesh, angle, estimable);
    exportAsPly(prefix + "_cv_l.ply", mesh, cvL, estimable);
    exportAsPly(prefix + "_ratio.ply", mesh, ratio, estimable);
    return 0;
}
