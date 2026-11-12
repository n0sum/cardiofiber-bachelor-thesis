// FiberView: a conduction field in format D, as two PLY files for MeshLab.
//
// Usage, from the repository root:
//   ./estimator/build/FiberView D.ply PREFIX
//   ./estimator/build/FiberView data/reference/uniform/truth.ply data/reference/uniform/truth
//
// PREFIX_ratio.ply   the mesh, quality per face = anisotropy ratio r = CV_L / CV_T = exp(2 alpha)
// PREFIX_fibers.ply  per face a thin triangle along the fiber, slightly above the face;
//                    its length is K * CV_L. View it in wireframe.
//
// Faces with a NaN in fx, fy, fz, mu or alpha have no value and are left out of both files.

#include "io.h"
#include <happly.h>
#include <glm/glm.hpp>
#include <cardiofiber/Mesh.h>

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

// Length of a fiber triangle per unit of CV_L (mm per mm/ms).
const double K = 0.6;
// Width of the base of a fiber triangle (mm).
const double WIDTH = 0.05;
// Height of a fiber triangle above its face, along the outward normal (mm).
const double HEIGHT = 0.05;

int main(int argc, char** argv) {
    if (argc != 3) {
        throw std::runtime_error("usage: FiberView D.ply PREFIX");
    }
    std::string inputFile = argv[1];
    std::string prefix = argv[2];

    TensorPly field = readTensorPly(inputFile);
    int numFaces = (int)field.faces.size();

    // --- which faces have a value: fiber, mu and alpha all finite
    std::vector<bool> hasValue(numFaces);
    for (int i = 0; i < numFaces; i++) {
        glm::dvec3 f = field.fibers[i];
        hasValue[i] = std::isfinite(f.x) && std::isfinite(f.y) && std::isfinite(f.z) &&
                      std::isfinite(field.mu[i]) && std::isfinite(field.alpha[i]);
    }

    // --- PREFIX_ratio.ply: anisotropy ratio per face
    std::vector<double> ratio(numFaces);
    for (int i = 0; i < numFaces; i++) {
        ratio[i] = std::exp(2.0 * field.alpha[i]);      // NaN where the face has no value
    }
    cardiofiber::Mesh mesh(field.vertices, field.faces, {});
    exportAsPly(prefix + "_ratio.ply", mesh, ratio, hasValue);

    // --- PREFIX_fibers.ply: one thin triangle per face with a value, three new vertices each
    std::vector<float> x, y, z;
    std::vector<std::vector<int>> triangles;
    for (int i = 0; i < numFaces; i++) {
        if (!hasValue[i]) continue;

        glm::dvec3 a = field.vertices[field.faces[i].a];
        glm::dvec3 b = field.vertices[field.faces[i].b];
        glm::dvec3 c = field.vertices[field.faces[i].c];
        glm::dvec3 center = (a + b + c) / 3.0;
        glm::dvec3 n = glm::normalize(glm::cross(b - a, c - a));   // outwards
        glm::dvec3 f = glm::normalize(field.fibers[i]);
        glm::dvec3 t = glm::cross(n, f);                            // across the fiber, in the face plane

        double cvL = std::exp(field.mu[i] + field.alpha[i]);
        double L = K * cvL;

        glm::dvec3 lifted = center + HEIGHT * n;
        glm::dvec3 corners[3] = {
            lifted - (L / 2.0) * f - (WIDTH / 2.0) * t,
            lifted - (L / 2.0) * f + (WIDTH / 2.0) * t,
            lifted + (L / 2.0) * f
        };

        int first = (int)x.size();
        for (int j = 0; j < 3; j++) {
            x.push_back(static_cast<float>(corners[j].x));
            y.push_back(static_cast<float>(corners[j].y));
            z.push_back(static_cast<float>(corners[j].z));
        }
        triangles.push_back({first, first + 1, first + 2});
    }

    happly::PLYData plyOut;
    plyOut.addElement("vertex", x.size());
    plyOut.getElement("vertex").addProperty<float>("x", x);
    plyOut.getElement("vertex").addProperty<float>("y", y);
    plyOut.getElement("vertex").addProperty<float>("z", z);
    plyOut.addElement("face", triangles.size());
    plyOut.getElement("face").addListProperty<int>("vertex_indices", triangles);
    plyOut.write(prefix + "_fibers.ply", happly::DataFormat::ASCII);
    return 0;
}
