#include <iostream>
#include <vector>
#include <glm/vec3.hpp>
#include <cardiofiber/Mesh.h>
#include "ExperimentConfig.h"
#include "io.h"
#include "GroundTruth.h"

// Usage, from the repository root:
//   ./estimator/build/Estimate experiments/a.conf [experiments/b.conf ...]
// One CSV of per-face errors per experiment, in sequence. The folder of output_csv must exist:
// experiments/subsets.py creates it when it writes the .conf.

int main(int argc, char** argv) {
    for (int i = 1; i < argc; ++i) {
        const ExperimentConfig cfg = loadExperimentConfig(argv[i]);
        std::cout << argv[i] << " -> " << cfg.outputCsv << "\n";
        const MultiLatPly in = readMultiLatPly(cfg.latPly, cfg.maps);
        const cardiofiber::Mesh m(in.vertices, in.faces, in.lats);
        const std::vector<std::vector<glm::dvec3>> gradients = m.computeGradients();
        const std::vector<cardiofiber::FaceAnisotropy> est = m.estimateAnisotropy(gradients, cfg.minSingularValue);
        const std::vector<FaceGroundTruth> gt = readGroundTruthPly(cfg.groundTruthPly);
        exportErrorsAsCsv(cfg.outputCsv, m, est, gt);
    }
    return 0;
}
