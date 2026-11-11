#pragma once
#include <string>
#include <vector>

// One experiment read from experiments/*.conf (format in experiments/README.md).
// Paths relative to the repository root: the executable is launched from there.
struct ExperimentConfig {
    std::string latPly, groundTruthPly, outputCsv;
    std::vector<int> maps;              // indices of the LAT maps (lat_K)
    double minSingularValue = 0.0;      
};

ExperimentConfig loadExperimentConfig(const std::string& path);
