#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>
#include <glm/vec3.hpp>
#include <cardiofiber/Mesh.h>
#include "ExperimentConfig.h"
#include "io.h"

// Usage, from the repository root:
//   ./estimator/build/Estimate experiments/a.conf [experiments/b.conf ...]
// For each experiment, in sequence: the LAT maps of lat_ply -> per-face estimate.
//   output_ply  the estimate in format D (mesh of lat_ply as it is; NaN on the non-estimable faces)
//   output_csv  face_id,cond, one row per face; cond empty where it is not computed
// The estimator never reads the ground truth. The folder of the outputs must exist:
// experiments/subsets.py creates it when it writes the .conf.

void writeCondCsv(const std::string& filename, const std::vector<cardiofiber::FaceAnisotropy>& est) {
    std::ofstream out(filename);
    if (!out) throw std::runtime_error("cannot open " + filename);
    out << std::scientific << std::setprecision(9);
    out << "face_id,cond\n";
    for (size_t i = 0; i < est.size(); ++i) {
        out << i << ',';
        if (est[i].cond > 0.0) out << est[i].cond;          // empty if not computed
        out << '\n';
    }
}

int main(int argc, char** argv) {
    for (int i = 1; i < argc; ++i) {
        const ExperimentConfig cfg = loadExperimentConfig(argv[i]);
        std::cout << argv[i] << " -> " << cfg.outputPly << "\n";
        const MultiLatPly in = readMultiLatPly(cfg.latPly, cfg.maps);
        const cardiofiber::Mesh m(in.vertices, in.faces, in.lats);
        const std::vector<std::vector<glm::dvec3>> gradients = m.computeGradients();
        const std::vector<cardiofiber::FaceAnisotropy> est = m.estimateAnisotropy(gradients, cfg.minSingularValue);

        const double nan = std::numeric_limits<double>::quiet_NaN();
        TensorPly out;
        out.vertices = in.vertices;
        out.faces = in.faces;
        for (const cardiofiber::FaceAnisotropy& a : est) {
            if (a.estimable) {
                out.fibers.push_back(a.fiber_dir);
                out.mu.push_back((std::log(a.CV_L) + std::log(a.CV_T)) / 2.0);
                out.alpha.push_back((std::log(a.CV_L) - std::log(a.CV_T)) / 2.0);
            } else {
                out.fibers.push_back(glm::dvec3(nan));
                out.mu.push_back(nan);
                out.alpha.push_back(nan);
            }
        }
        writeTensorPly(cfg.outputPly, out);
        writeCondCsv(cfg.outputCsv, est);
    }
    return 0;
}
