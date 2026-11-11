#pragma once
#include <string>
#include <vector>
#include <glm/vec3.hpp>
#include <cardiofiber/Mesh.h>

struct FaceGroundTruth {
    glm::dvec3 fiber_dir = glm::dvec3(0);
    double CV_L = 0.0;              
    double CV_T = 0.0;              
};

std::vector<FaceGroundTruth> readGroundTruthPly(const std::string& filename);

void exportErrorsAsCsv(const std::string& filename,
                       const cardiofiber::Mesh& mesh, const std::vector<cardiofiber::FaceAnisotropy>& est,
                       const std::vector<FaceGroundTruth>& gt);
