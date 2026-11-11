#pragma once
#include <string>
#include <vector>
#include <glm/vec3.hpp>
#include <cardiofiber/Mesh.h>

struct MultiLatPly {
    std::vector<glm::dvec3> vertices;
    std::vector<cardiofiber::Face> faces;
    std::vector<std::vector<double>> lats;      //[map][vertex]
};

MultiLatPly readMultiLatPly(const std::string& filename, const std::vector<int>& maps);

void exportAsPly(const std::string& filename, const cardiofiber::Mesh& mesh,
                 const std::vector<double>& faceQuality, const std::vector<bool>& writeFace);
