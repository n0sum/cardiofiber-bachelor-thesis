// Simulate: activation maps (LAT) from a conduction field in format D.
//
// Solves the anisotropic eikonal equation on the triangle mesh with the FIM
// discretization of Fu et al. 2011, in its simple form: no active list, sweeps
// over all faces until nothing changes. One map per row of SITES.csv.
// With SIGMA > 0, Gaussian noise is added to the LATs before writing: to every
// vertex of every map a number drawn from a normal distribution with mean 0 and
// standard deviation SIGMA (ms), independent per vertex and per map, also on the
// vertices of the stimulus sphere. One generator with the fixed seed NOISE_SEED,
// so the same SIGMA gives the same file. With SIGMA = 0 nothing is drawn.
//
// Usage, from the repository root:
//   ./estimator/build/Simulate TRUTH.ply SITES.csv LAT.ply SIGMA
//   ./estimator/build/Simulate data/reference/uniform/truth.ply simulator/stage1/pacing_sites.csv data/reference/uniform/lat.ply 0
//
// TRUTH.ply: format D (per vertex x y z in mm; per face vertex_indices,
//            fiber fx fy fz, mu, alpha; CV_L = exp(mu + alpha), CV_T = exp(mu - alpha) in mm/ms).
//            The faces are used as they are, nothing is flipped here.
// SITES.csv: header cell_id, one face index per row; row k is map k.
// LAT.ply:   per vertex x y z (mm) and lat_0 ... lat_K (ms), faces vertex_indices.
// SIGMA:     standard deviation of the noise on the LATs, ms; 0 = no noise.

#include <happly.h>
#include <glm/glm.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

// Vertices within this distance (mm) of the barycenter of the pacing face start at LAT 0.
const double STIMULUS_RADIUS = 3.0;

// Seed of the generator of the noise on the LATs: fixed, so the runs are reproducible.
const unsigned int NOISE_SEED = 12345;

const double INF = std::numeric_limits<double>::infinity();

// Reads the face indices of SITES.csv: a header line, then one integer per line.
std::vector<int> readSites(const std::string& filename) {
    std::ifstream in(filename);
    if (!in) throw std::runtime_error("cannot open " + filename);
    std::vector<int> sites;
    std::string line;
    std::getline(in, line);                     // header: cell_id
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        sites.push_back(std::stoi(line));
    }
    return sites;
}

// Travel time of the segment d inside a face with metric W: sqrt(d^T W d).
double travelTime(const glm::dmat3& W, const glm::dvec3& d) {
    return std::sqrt(glm::dot(d, W * d));
}

// Smallest arrival time at vertex C from the other two vertices A and B of a face.
double updateFromFace(const glm::dvec3& A, const glm::dvec3& B, const glm::dvec3& C,
                      double TA, double TB, const glm::dmat3& W) {
    glm::dvec3 u = C - A;
    glm::dvec3 v = C - B;

    // from the edges
    double best = INF;
    if (std::isfinite(TA)) best = std::min(best, TA + travelTime(W, u));
    if (std::isfinite(TB)) best = std::min(best, TB + travelTime(W, v));

    // from the triangle: the front reaches C from P = A + s*e on edge AB,
    // with the time linear along AB
    if (std::isfinite(TA) && std::isfinite(TB)) {
        glm::dvec3 e = B - A;
        double delta = TB - TA;
        double h1 = glm::dot(e, W * e);
        double h2 = glm::dot(e, W * u);
        double h3 = glm::dot(u, W * u);
        if (delta * delta < h1) {
            double s = (h2 - delta * std::sqrt((h1 * h3 - h2 * h2) / (h1 - delta * delta))) / h1;
            if (s >= 0.0 && s <= 1.0) {
                double T = TA + s * delta + std::sqrt(h1 * s * s - 2.0 * h2 * s + h3);
                best = std::min(best, T);
            }
        }
    }
    return best;
}

int main(int argc, char** argv) {
    if (argc != 5) {
        throw std::runtime_error("usage: Simulate TRUTH.ply SITES.csv LAT.ply SIGMA");
    }
    std::string truthFile = argv[1];
    std::string sitesFile = argv[2];
    std::string latFile = argv[3];
    double sigma = std::stod(argv[4]);

    // --- read the conduction field
    happly::PLYData plyIn(truthFile);
    happly::Element& vertexElement = plyIn.getElement("vertex");
    happly::Element& faceElement = plyIn.getElement("face");
    std::vector<double> x = vertexElement.getProperty<double>("x");
    std::vector<double> y = vertexElement.getProperty<double>("y");
    std::vector<double> z = vertexElement.getProperty<double>("z");
    std::vector<std::vector<int>> faces = faceElement.getListProperty<int>("vertex_indices");
    std::vector<double> fx = faceElement.getProperty<double>("fx");
    std::vector<double> fy = faceElement.getProperty<double>("fy");
    std::vector<double> fz = faceElement.getProperty<double>("fz");
    std::vector<double> mu = faceElement.getProperty<double>("mu");
    std::vector<double> alpha = faceElement.getProperty<double>("alpha");

    int numVertices = (int)x.size();
    int numFaces = (int)faces.size();

    std::vector<glm::dvec3> vertices(numVertices);
    for (int i = 0; i < numVertices; i++) {
        vertices[i] = glm::dvec3(x[i], y[i], z[i]);
    }

    // --- metric per face: W = f f^T / CV_L^2 + (I - f f^T) / CV_T^2
    std::vector<glm::dmat3> metric(numFaces);
    for (int i = 0; i < numFaces; i++) {
        double cvL = std::exp(mu[i] + alpha[i]);
        double cvT = std::exp(mu[i] - alpha[i]);
        glm::dvec3 f(fx[i], fy[i], fz[i]);
        glm::dmat3 ffT = glm::outerProduct(f, f);
        glm::dmat3 identity(1.0);
        metric[i] = ffT / (cvL * cvL) + (identity - ffT) / (cvT * cvT);
    }

    std::vector<int> sites = readSites(sitesFile);
    std::vector<std::vector<double>> lats;

    // --- generator of the noise, one for all the maps
    std::mt19937 generator(NOISE_SEED);
    std::normal_distribution<double> noise(0.0, 1.0);

    for (int k = 0; k < (int)sites.size(); k++) {
        // --- start: LAT 0 within STIMULUS_RADIUS of the barycenter of the pacing face
        const std::vector<int>& site = faces[sites[k]];
        glm::dvec3 center = (vertices[site[0]] + vertices[site[1]] + vertices[site[2]]) / 3.0;
        std::vector<double> T(numVertices);
        for (int i = 0; i < numVertices; i++) {
            if (glm::length(vertices[i] - center) <= STIMULUS_RADIUS) T[i] = 0.0;
            else T[i] = INF;
        }

        // --- sweeps over all faces, updating in place, until a sweep lowers nothing
        bool changed = true;
        while (changed) {
            changed = false;
            for (int i = 0; i < numFaces; i++) {
                const std::vector<int>& face = faces[i];
                for (int j = 0; j < 3; j++) {
                    int c = face[j];
                    int a = face[(j + 1) % 3];
                    int b = face[(j + 2) % 3];
                    double candidate = updateFromFace(vertices[a], vertices[b], vertices[c],
                                                      T[a], T[b], metric[i]);
                    if (candidate < T[c]) {
                        T[c] = candidate;
                        changed = true;
                    }
                }
            }
        }

        // --- check: the mesh is one component, every vertex must be reached
        int unreached = 0;
        for (int i = 0; i < numVertices; i++) {
            if (std::isinf(T[i])) unreached++;
        }
        if (unreached > 0) {
            throw std::runtime_error("map " + std::to_string(k) + ": " + std::to_string(unreached) +
                                     " vertices not reached");
        }

        // --- noise: to every vertex a number from a Gaussian with mean 0 and standard deviation sigma
        if (sigma > 0.0) {
            for (int i = 0; i < numVertices; i++) {
                T[i] = T[i] + sigma * noise(generator);
            }
        }
        lats.push_back(T);
    }

    // --- write the maps
    happly::PLYData plyOut;
    plyOut.addElement("vertex", numVertices);
    plyOut.getElement("vertex").addProperty<double>("x", x);
    plyOut.getElement("vertex").addProperty<double>("y", y);
    plyOut.getElement("vertex").addProperty<double>("z", z);
    for (int k = 0; k < (int)lats.size(); k++) {
        plyOut.getElement("vertex").addProperty<double>("lat_" + std::to_string(k), lats[k]);
    }
    plyOut.addElement("face", numFaces);
    plyOut.getElement("face").addListProperty<int>("vertex_indices", faces);
    plyOut.write(latFile, happly::DataFormat::Binary);
    return 0;
}
