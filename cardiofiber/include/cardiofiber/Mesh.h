//
// Created by gian on 12/19/25.
//

#pragma once
#include <vector>
#include <glm/vec3.hpp>

namespace cardiofiber {

struct Face { int a, b, c; };

// Anisotropy estimate of one face.
// Local tangent basis, deterministic (depends only on the stored vertex order):
//   e1 = p_b - p_a,  e2 = p_c - p_a,  n = normalize(cross(e1, e2)),  u = normalize(e1),  w = cross(n, u)
// The tensor (mm^2/ms^2) is estimated in this basis; fiber_dir is mapped back to 3D.
struct FaceAnisotropy {
    double CV_L = 0.0;              
    double CV_T = 0.0;              
    glm::dvec3 fiber_dir = glm::dvec3(0); 
    double cond = 0.0;              
    bool estimable = false;      
};

class Mesh{
    std::vector<glm::dvec3> vPos;
    std::vector<Face> fIndices;
    std::vector<std::vector<double>> vQualities;   //[map][vertex]

public:
    Mesh(const std::vector<glm::dvec3>& vertices, const std::vector<Face>& faces,
         const std::vector<std::vector<double>>& lats);

    const std::vector<glm::dvec3>& vertices() const;
    const std::vector<Face>& faces() const;

    std::vector<glm::dvec3> areaVectors() const;
    // Per-vertex normals, area-weighted
    std::vector<glm::dvec3> computeVertexNormals() const;

    // gradients[map][face].
    std::vector<std::vector<glm::dvec3>> computeGradients() const;

    // The gradients can come from computeGradients() or from outside: in 3D, tangent to the face
    FaceAnisotropy estimateFace(int i, const std::vector<std::vector<glm::dvec3>>& gradients,
                                double minSingularValue) const;

    std::vector<FaceAnisotropy> estimateAnisotropy(const std::vector<std::vector<glm::dvec3>>& gradients,
                                                   double minSingularValue) const;

private:
    std::vector<glm::dvec3> distributeNormalOnVertex(const std::vector<glm::dvec3> &fNorm) const;
    glm::dvec3 faceGradientOf(int faceIndex, int mapIndex) const;
};

}
