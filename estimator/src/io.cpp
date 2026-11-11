#include "io.h"
#include <happly.h>
#include <string>
#include <vector>
#include <stdexcept>
#include <cmath>

MultiLatPly readMultiLatPly(const std::string& filename, const std::vector<int>& maps) {
    MultiLatPly in;

    happly::PLYData plyIn(filename);
    auto& v = plyIn.getElement("vertex");

    // float32 in the file; happly converts to double, exactly
    std::vector<double> x = v.getProperty<double>("x");
    std::vector<double> y = v.getProperty<double>("y");
    std::vector<double> z = v.getProperty<double>("z");

    in.vertices.resize(x.size());
    for (size_t i = 0; i < x.size(); ++i) in.vertices[i] = glm::dvec3(x[i], y[i], z[i]);

    const std::vector<std::vector<int>> prov_faces = plyIn.getElement("face").getListProperty<int>("vertex_indices");
    in.faces.reserve(prov_faces.size());
    for (const std::vector<int>& f : prov_faces) in.faces.emplace_back(f[0], f[1], f[2]);

    for (int k : maps) in.lats.push_back(v.getProperty<double>("lat_" + std::to_string(k)));
    return in;
}

void exportAsPly(const std::string &filename, const cardiofiber::Mesh& mesh, const std::vector<double>& faceQuality,
                 const std::vector<bool>& writeFace) {
    const std::vector<glm::dvec3>& vPos = mesh.vertices();
    const std::vector<cardiofiber::Face>& faces = mesh.faces();
    if (faceQuality.size() != faces.size())
        throw std::invalid_argument("exportAsPly: " + std::to_string(faceQuality.size()) + " quality values for " +
                                    std::to_string(faces.size()) + " faces");
    if (writeFace.size() != faces.size())
        throw std::invalid_argument("exportAsPly: " + std::to_string(writeFace.size()) + " mask values for " +
                                    std::to_string(faces.size()) + " faces");
    const std::vector<glm::dvec3> vNorm = mesh.computeVertexNormals();
    // the file is float32: conversion from double only here
    std::vector<float> x(vPos.size()),y(vPos.size()),z(vPos.size()),nx(vPos.size()),ny(vPos.size()),nz(vPos.size());
    for (int i = 0; i < vPos.size(); i++) {
        x[i] = static_cast<float>(vPos[i].x);
        y[i] = static_cast<float>(vPos[i].y);
        z[i] = static_cast<float>(vPos[i].z);
        nx[i] = static_cast<float>(vNorm[i].x);
        ny[i] = static_cast<float>(vNorm[i].y);
        nz[i] = static_cast<float>(vNorm[i].z);
    }
    happly::PLYData plyOut;
    plyOut.addElement("vertex", vPos.size());
    plyOut.getElement("vertex").addProperty<float>("x", x);
    plyOut.getElement("vertex").addProperty<float>("y", y);
    plyOut.getElement("vertex").addProperty<float>("z", z);
    plyOut.getElement("vertex").addProperty<float>("nx", nx);
    plyOut.getElement("vertex").addProperty<float>("ny", ny);
    plyOut.getElement("vertex").addProperty<float>("nz", nz);
    std::vector<std::vector<int>> vertex_indices;
    std::vector<float> quality;
    for (size_t i = 0; i < faces.size(); ++i) {
        if (!writeFace[i]) continue;
        if (!std::isfinite(faceQuality[i]))
            throw std::runtime_error(filename + ": non-finite quality on face_id " + std::to_string(i));
        vertex_indices.push_back({faces[i].a, faces[i].b, faces[i].c});
        quality.push_back(static_cast<float>(faceQuality[i]));
    }
    plyOut.addElement("face", vertex_indices.size());
    plyOut.getElement("face").addListProperty<int>("vertex_indices", vertex_indices);
    plyOut.getElement("face").addProperty<float>("quality", quality);
    plyOut.write(filename, happly::DataFormat::ASCII);
}
