// ErrorMaps: the log-Euclidean error of an estimate as a PLY for MeshLab, one RGB color per face.
//
// Usage, from the repository root:
//   ./estimator/build/ErrorMaps MESH.ply NUMBERS.txt OUT.ply
// MESH.ply: a PLY of the mesh (vertices and faces; e.g. a lat.ply), faces used as they are, nothing is flipped.
// NUMBERS.txt: one row per face, "estimable error" (experiments/error_maps.py), the log-Euclidean
// distance between estimated and true tensor; "0 0" on the non-estimable faces.
// Writes OUT.ply (writeColoredPly, io.h): all the vertices and all the faces of MESH.ply, in order; the
// non-estimable faces gray, the others colored by the error, green (0) -> yellow (E_MAX / 2) -> red
// (E_MAX and above). In MeshLab: Render -> Color -> Per Face.

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "io.h"

// Color of the non-estimable faces.
const std::array<unsigned char, 3> GRAY = {160, 160, 160};
// Error at which the color saturates to red: 2 sqrt 2 ln 2 ~ 1.96, the error of a fiber rotated by
// 90 degrees with r = 2 (speeds exact).
const double E_MAX = 2.0 * std::sqrt(2.0) * std::log(2.0);

// green -> yellow -> red: t = min(E / E_MAX, 1), red = 255 min(1, 2t), green = 255 min(1, 2(1 - t)), blue = 0.
std::array<unsigned char, 3> errorColor(double error) {
    const double t = std::min(error / E_MAX, 1.0);
    return {static_cast<unsigned char>(std::lround(255.0 * std::min(1.0, 2.0 * t))),
            static_cast<unsigned char>(std::lround(255.0 * std::min(1.0, 2.0 * (1.0 - t)))),
            0};
}

int main(int argc, char** argv) {
    if (argc != 4) {
        throw std::runtime_error("usage: ErrorMaps MESH.ply NUMBERS.txt OUT.ply");
    }
    const MultiLatPly in = readMultiLatPly(argv[1], {});

    std::vector<std::array<unsigned char, 3>> colors;
    std::ifstream numbers(argv[2]);
    int e;
    double v;
    while (numbers >> e >> v) {
        if (e == 1 && !std::isfinite(v))
            throw std::runtime_error(std::string(argv[2]) + ": non-finite error on face_id " +
                                     std::to_string(colors.size()));
        colors.push_back(e == 1 ? errorColor(v) : GRAY);
    }

    writeColoredPly(argv[3], in.vertices, in.faces, colors);
    return 0;
}
