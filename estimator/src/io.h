#pragma once
#include <array>
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

// Conduction field in format D: CV_L = exp(mu + alpha), CV_T = exp(mu - alpha), mm/ms.
// Faces as in the file (normals already outwards); NaN on faces without a value.
struct TensorPly {
    std::vector<glm::dvec3> vertices;
    std::vector<cardiofiber::Face> faces;
    std::vector<glm::dvec3> fibers;   // per face
    std::vector<double> mu;           // per face
    std::vector<double> alpha;        // per face
};

TensorPly readTensorPly(const std::string& filename);
// Text PLY (format ascii 1.0), with the comment line of format D; NaN values are written as they are ("nan").
void writeTensorPly(const std::string& filename, const TensorPly& out);

// Text PLY (format ascii 1.0) for MeshLab: per vertex "property float x/y/z" (no normals), per face the
// triangle and its color "property uchar red/green/blue". Faces written as given, in order; nothing is flipped.
void writeColoredPly(const std::string& filename, const std::vector<glm::dvec3>& vertices,
                     const std::vector<cardiofiber::Face>& faces,
                     const std::vector<std::array<unsigned char, 3>>& colors);
