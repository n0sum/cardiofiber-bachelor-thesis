//
// Created by gian on 12/19/25.
//

#include "cardiofiber/Mesh.h"
#include <vector>
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <glm/geometric.hpp>
#include <glm/vec2.hpp>
#include <glm/mat2x2.hpp>

using namespace glm;

namespace cardiofiber {


std::vector<dvec3> normalizeVectors(const std::vector<dvec3> &vectors) {
    std::vector<dvec3> result(vectors.size());
    for (int i = 0; i<vectors.size(); ++i) {
        result[i] = normalize(vectors[i]);
    }
    return result;
}

// Structural thresholds of the anisotropy estimate. The threshold on sigma_min
// comes from the experiment file.
constexpr int    kMinGradients       = 3;       // below: underdetermined system
// Numerical zero of lambda2 = CV_T^2, relative to lambda1 = CV_L^2: a test of numerical zero, not a
// physical threshold. It means CV_T / CV_L > 1e-4. Without a tolerance the exactly flat ellipses
// (tensor of rank 1, lambda2 = 0 in exact arithmetic) pass or not depending on the rounding.
constexpr double kMinEigenRatio      = 1e-8;

struct LsqResult {
    std::array<double,3> m{};
    double sigmaMax = 0.0, sigmaMin = 0.0;
};

// Least squares R m = b with R (K x 3) via one-sided Jacobi SVD: rotates the columns of R
// until they are orthogonal, so that R V = U S with S = column norms.
// m = V S^-1 U^T b = sum_j v_j (u_j . b) / sigma_j. m is always computed: with a zero sigma it
// comes out inf or NaN, so the caller checks sigmaMin before using it.
LsqResult svdLeastSquares(std::vector<std::array<double,3>> R, const std::vector<double> &b) {
    LsqResult out;
    const int K = R.size();
    double V[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
    for (int sweep = 0; sweep < 60; ++sweep) {
        bool rotated = false;
        for (int p = 0; p < 2; ++p) for (int q = p + 1; q < 3; ++q) {
            double alpha = 0, beta = 0, gamma = 0;
            for (int i = 0; i < K; ++i) {
                alpha += R[i][p]*R[i][p];
                beta  += R[i][q]*R[i][q];
                gamma += R[i][p]*R[i][q];
            }
            if (std::abs(gamma) <= 1e-15 * std::sqrt(alpha*beta)) continue;
            rotated = true;
            const double zeta = (beta - alpha) / (2.0*gamma);
            const double t = (zeta >= 0 ? 1.0 : -1.0) / (std::abs(zeta) + std::sqrt(1.0 + zeta*zeta));
            const double c = 1.0 / std::sqrt(1.0 + t*t), s = c*t;
            for (int i = 0; i < K; ++i) {
                const double rp = R[i][p], rq = R[i][q];
                R[i][p] = c*rp - s*rq;
                R[i][q] = s*rp + c*rq;
            }
            for (int i = 0; i < 3; ++i) {
                const double vp = V[i][p], vq = V[i][q];
                V[i][p] = c*vp - s*vq;
                V[i][q] = s*vp + c*vq;
            }
        }
        if (!rotated) break;
    }
    double sigma[3];
    for (int j = 0; j < 3; ++j) {
        double n2 = 0;
        for (int i = 0; i < K; ++i) n2 += R[i][j]*R[i][j];
        sigma[j] = std::sqrt(n2);
    }
    out.sigmaMax = std::max({sigma[0], sigma[1], sigma[2]});
    out.sigmaMin = std::min({sigma[0], sigma[1], sigma[2]});
    out.m = {0,0,0};
    for (int j = 0; j < 3; ++j) {
        double ub = 0;                                   // u_j . b, with u_j = column_j / sigma_j
        for (int i = 0; i < K; ++i) ub += R[i][j]*b[i];
        ub /= sigma[j];
        for (int r = 0; r < 3; ++r) out.m[r] += V[r][j]*ub/sigma[j];
    }
    return out;
}

// processing methods

Mesh::Mesh(const std::vector<dvec3>& vertices, const std::vector<Face>& faces,
           const std::vector<std::vector<double>>& lats)
    : vPos(vertices), fIndices(faces), vQualities(lats) {
    for (size_t k = 0; k < vQualities.size(); ++k)
        if (vQualities[k].size() != vPos.size())
            throw std::invalid_argument("Mesh: map " + std::to_string(k) + " has " +
                                        std::to_string(vQualities[k].size()) + " LATs for " +
                                        std::to_string(vPos.size()) + " vertices");
}

const std::vector<dvec3>& Mesh::vertices() const { return vPos; }

const std::vector<Face>& Mesh::faces() const { return fIndices; }

std::vector<dvec3> Mesh::areaVectors() const {
    std::vector<dvec3>  result(fIndices.size() );
   for (int i = 0; Face face : fIndices) {
        dvec3 e1 = vPos[face.b] - vPos[face.a];
        dvec3 e2 = vPos[face.c] - vPos[face.a];
        result[i++] = cross(e1,e2);
    }
    return result;
}

std::vector<dvec3> Mesh::distributeNormalOnVertex(const std::vector<dvec3> &fNorm) const {
    std::vector<dvec3> result(vPos.size());
    for (int i = 0; i < fIndices.size(); ++i) {
        result[fIndices[i].a] += fNorm[i];
        result[fIndices[i].b] += fNorm[i];
        result[fIndices[i].c] += fNorm[i];
    }
    return result;
}

dvec3 Mesh::faceGradientOf(int faceIndex, int mapIndex) const {
    Face face = fIndices[faceIndex];
    dvec3 e1 = vPos[face.b] - vPos[face.a];
    dvec3 e2 = vPos[face.c] - vPos[face.a];
    dmat2 A(
        dot(e1,e1),dot(e1,e2),
        dot(e1,e2), dot(e2,e2)
        );
    dvec2 k(vQualities[mapIndex][face.b]-vQualities[mapIndex][face.a], vQualities[mapIndex][face.c]-vQualities[mapIndex][face.a]);
    dmat2 inv = inverse(A);
    dvec2 solution = inv*k;
    return (solution.x * e1 + solution.y * e2);
}

std::vector<dvec3> Mesh::computeVertexNormals() const {
    return normalizeVectors(distributeNormalOnVertex(areaVectors()));
}

std::vector<std::vector<dvec3>> Mesh::computeGradients() const {
    std::vector<std::vector<dvec3>> gradients(vQualities.size(), std::vector<dvec3>(fIndices.size()));
    for (int k = 0; k < vQualities.size(); ++k)
        for (int i = 0; i < fIndices.size(); ++i)
            gradients[k][i] = faceGradientOf(i, k);
    return gradients;
}

FaceAnisotropy Mesh::estimateFace(int i, const std::vector<std::vector<dvec3>>& gradients,
                                  double minSingularValue) const {
    FaceAnisotropy fa;
    const int K = gradients.size();
    const Face &face = fIndices[i];

    if (K < kMinGradients) return fa;

    const dvec3 e1 = vPos[face.b] - vPos[face.a];
    const dvec3 e2 = vPos[face.c] - vPos[face.a];
    const dvec3 a2 = cross(e1, e2);
    const dvec3 n = a2 / length(a2);
    const dvec3 u = e1 / length(e1);
    const dvec3 w = cross(n, u);

    std::vector<std::array<double,3>> R(K);
    const std::vector<double> b(K, 1.0);
    for (int k = 0; k < K; ++k) {
        const dvec3 g = gradients[k][i];
        const double g1 = dot(g, u), g2 = dot(g, w);
        R[k] = {g1*g1, std::sqrt(2.0)*g1*g2, g2*g2};
    }

    const LsqResult lsq = svdLeastSquares(R, b);
    // the form !(x > threshold) also discards NaN
    if (!(lsq.sigmaMin > minSingularValue)) return fa;
    fa.cond = lsq.sigmaMax / lsq.sigmaMin;

    const double M11 = lsq.m[0], M12 = lsq.m[1] / std::sqrt(2.0), M22 = lsq.m[2];
    const double sum = M11 + M22;
    const double diff = std::sqrt((M11 - M22)*(M11 - M22) + (2.0*M12)*(2.0*M12));
    const double lambda1 = (sum + diff) / 2.0;           // CV_L^2
    const double lambda2 = (sum - diff) / 2.0;           // CV_T^2
    // CV_T^2 must be positive, up to the numerical zero kMinEigenRatio; 
    if (!(lambda2 > kMinEigenRatio * lambda1)) return fa;

    const double theta = 0.5 * std::atan2(2.0*M12, M11 - M22);
    fa.CV_L = std::sqrt(lambda1);
    fa.CV_T = std::sqrt(lambda2);
    fa.fiber_dir = std::cos(theta)*u + std::sin(theta)*w;
    fa.estimable = true;
    return fa;
}

std::vector<FaceAnisotropy> Mesh::estimateAnisotropy(const std::vector<std::vector<dvec3>>& gradients,
                                                     double minSingularValue) const {
    for (size_t k = 0; k < gradients.size(); ++k)
        if (gradients[k].size() != fIndices.size())
            throw std::invalid_argument("estimateAnisotropy: map " + std::to_string(k) + " has " +
                                        std::to_string(gradients[k].size()) + " gradients for " +
                                        std::to_string(fIndices.size()) + " faces");
    std::vector<FaceAnisotropy> est(fIndices.size());
    for (int i = 0; i < fIndices.size(); ++i) est[i] = estimateFace(i, gradients, minSingularValue);
    return est;
}

}
