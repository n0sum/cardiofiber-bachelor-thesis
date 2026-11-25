// FiberView: a conduction field in format D as one PLY for MeshLab, the fibers drawn as ellipses.
//
// Usage, from the repository root:
//   ./estimator/build/FiberView D.ply OUT.ply
//   ./estimator/build/FiberView data/reference/painted/truth.ply output/mesh/ellipses/truth_painted_ellipses.ply
//
// OUT.ply (writeColoredPly, io.h), one RGB color per face, in the same file:
// - the mesh of D.ply, all its faces gray, in order (nothing is flipped);
// - for each face with a value (fx, fy, fz, mu, alpha all finite) an ellipse, a fan of N_SEGMENTS
//   triangles H above the face along its normal: semi-axis RHO along the fiber and RHO / r across it,
//   r = CV_L / CV_T = exp(2 alpha). The long axis is fixed: every ellipse lies inside the circle of
//   radius RHO, a circle with alpha = 0, a needle 2 RHO long with large alpha. Colored by the mean
//   speed v = exp(mu): blue (V_MIN and below) -> red (V_MAX and above).
// Faces without a value: only gray, no ellipse. In MeshLab: Render -> Color -> Per Face.

#include "io.h"

#include <glm/glm.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

// Color of the mesh faces.
const std::array<unsigned char, 3> GRAY = {160, 160, 160};
// Semi-axis of an ellipse along the fiber (mm); across it RHO / r.
const double RHO = 0.2;
// Height of an ellipse above its face, along the outward normal (mm).
const double H = 0.05;
// Triangles per ellipse.
const int N_SEGMENTS = 16;
// Mean speed v = exp(mu) (mm/ms) mapped to blue and to red; saturated outside.
const double V_MIN = 0.5;
const double V_MAX = 1.0;

// blue -> red: s = clamp((v - V_MIN) / (V_MAX - V_MIN), 0, 1), red = 255 s, green = 0, blue = 255 (1 - s).
std::array<unsigned char, 3> speedColor(double v) {
    const double s = std::clamp((v - V_MIN) / (V_MAX - V_MIN), 0.0, 1.0);
    return {static_cast<unsigned char>(std::lround(255.0 * s)), 0,
            static_cast<unsigned char>(std::lround(255.0 * (1.0 - s)))};
}

int main(int argc, char** argv) {
    if (argc != 3) {
        throw std::runtime_error("usage: FiberView D.ply OUT.ply");
    }
    const TensorPly field = readTensorPly(argv[1]);
    const int numFaces = (int)field.faces.size();

    // the mesh, all gray
    std::vector<glm::dvec3> vertices = field.vertices;
    std::vector<cardiofiber::Face> faces = field.faces;
    std::vector<std::array<unsigned char, 3>> colors(numFaces, GRAY);

    // one ellipse per face with a value: center + N_SEGMENTS rim vertices, N_SEGMENTS triangles
    for (int i = 0; i < numFaces; i++) {
        const glm::dvec3 fiber = field.fibers[i];
        if (!(std::isfinite(fiber.x) && std::isfinite(fiber.y) && std::isfinite(fiber.z) &&
              std::isfinite(field.mu[i]) && std::isfinite(field.alpha[i])))
            continue;

        const glm::dvec3 a = field.vertices[field.faces[i].a];
        const glm::dvec3 b = field.vertices[field.faces[i].b];
        const glm::dvec3 c = field.vertices[field.faces[i].c];
        const glm::dvec3 n = glm::normalize(glm::cross(b - a, c - a));      // outwards
        const glm::dvec3 f = glm::normalize(fiber - glm::dot(fiber, n) * n); // fiber in the face plane
        const glm::dvec3 t = glm::cross(n, f);                                // across the fiber
        const glm::dvec3 center = (a + b + c) / 3.0 + H * n;
        const double minor = RHO / std::exp(2.0 * field.alpha[i]);            // RHO / r
        const std::array<unsigned char, 3> color = speedColor(std::exp(field.mu[i]));

        const int first = (int)vertices.size();
        vertices.push_back(center);
        for (int j = 0; j < N_SEGMENTS; j++) {
            const double theta = 2.0 * std::numbers::pi * j / N_SEGMENTS;
            vertices.push_back(center + RHO * std::cos(theta) * f + minor * std::sin(theta) * t);
        }
        // (center, j, j + 1): normal along cross(f, t) = n, outwards
        for (int j = 0; j < N_SEGMENTS; j++) {
            faces.push_back({first, first + 1 + j, first + 1 + (j + 1) % N_SEGMENTS});
            colors.push_back(color);
        }
    }

    writeColoredPly(argv[2], vertices, faces, colors);
    return 0;
}
