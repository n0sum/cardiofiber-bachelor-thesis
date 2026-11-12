#include "ExperimentConfig.h"
#include <fstream>
#include <map>
#include <sstream>

std::vector<int> parseMapList(const std::string& v) {
    std::vector<int> result;
    std::stringstream ss(v);
    for (std::string tok; std::getline(ss, tok, ',');) result.push_back(std::stoi(tok));
    return result;
}

ExperimentConfig loadExperimentConfig(const std::string& path) {
    std::map<std::string, std::string> kv;
    std::ifstream in(path);
    for (std::string line; std::getline(in, line);) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        std::string key, equals, value;
        ss >> key >> equals >> value;
        kv[key] = value;
    }

    ExperimentConfig c;
    c.latPly = kv.at("lat_ply");
    c.maps = parseMapList(kv.at("maps"));
    c.minSingularValue = std::stod(kv.at("mask.min_singular_value"));
    c.outputPly = kv.at("output_ply");
    c.outputCsv = kv.at("output_csv");
    return c;
}
