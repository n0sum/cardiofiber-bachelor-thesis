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

    // double in lat.ply (Simulate)
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

TensorPly readTensorPly(const std::string& filename) {
    TensorPly in;

    happly::PLYData plyIn(filename);
    auto& v = plyIn.getElement("vertex");
    auto& f = plyIn.getElement("face");

    std::vector<double> x = v.getProperty<double>("x");
    std::vector<double> y = v.getProperty<double>("y");
    std::vector<double> z = v.getProperty<double>("z");

    in.vertices.resize(x.size());
    for (size_t i = 0; i < x.size(); ++i) in.vertices[i] = glm::dvec3(x[i], y[i], z[i]);

    const std::vector<std::vector<int>> prov_faces = f.getListProperty<int>("vertex_indices");
    in.faces.reserve(prov_faces.size());
    for (const std::vector<int>& face : prov_faces) in.faces.emplace_back(face[0], face[1], face[2]);

    std::vector<double> fx = f.getProperty<double>("fx");
    std::vector<double> fy = f.getProperty<double>("fy");
    std::vector<double> fz = f.getProperty<double>("fz");
    in.fibers.resize(fx.size());
    for (size_t i = 0; i < fx.size(); ++i) in.fibers[i] = glm::dvec3(fx[i], fy[i], fz[i]);

    in.mu = f.getProperty<double>("mu");
    in.alpha = f.getProperty<double>("alpha");
    return in;
}

void writeTensorPly(const std::string& filename, const TensorPly& out) {
    std::vector<double> x(out.vertices.size()), y(out.vertices.size()), z(out.vertices.size());
    for (size_t i = 0; i < out.vertices.size(); ++i) {
        x[i] = out.vertices[i].x;
        y[i] = out.vertices[i].y;
        z[i] = out.vertices[i].z;
    }
    std::vector<std::vector<int>> vertex_indices(out.faces.size());
    std::vector<double> fx(out.faces.size()), fy(out.faces.size()), fz(out.faces.size());
    for (size_t i = 0; i < out.faces.size(); ++i) {
        vertex_indices[i] = {out.faces[i].a, out.faces[i].b, out.faces[i].c};
        fx[i] = out.fibers[i].x;
        fy[i] = out.fibers[i].y;
        fz[i] = out.fibers[i].z;
    }

    happly::PLYData plyOut;
    plyOut.comments.push_back("mu = (ln CV_L + ln CV_T) / 2, alpha = (ln CV_L - ln CV_T) / 2, CV in mm/ms, natural log");
    plyOut.addElement("vertex", out.vertices.size());
    auto& v = plyOut.getElement("vertex");
    v.addProperty<double>("x", x);
    v.addProperty<double>("y", y);
    v.addProperty<double>("z", z);
    plyOut.addElement("face", out.faces.size());
    auto& f = plyOut.getElement("face");
    f.addListProperty<int>("vertex_indices", vertex_indices);
    f.addProperty<double>("fx", fx);
    f.addProperty<double>("fy", fy);
    f.addProperty<double>("fz", fz);
    f.addProperty<double>("mu", out.mu);
    f.addProperty<double>("alpha", out.alpha);
    plyOut.write(filename, happly::DataFormat::ASCII);
}

void writeColoredPly(const std::string& filename, const std::vector<glm::dvec3>& vertices,
                     const std::vector<cardiofiber::Face>& faces,
                     const std::vector<std::array<unsigned char, 3>>& colors) {
    if (colors.size() != faces.size())
        throw std::invalid_argument("writeColoredPly: " + std::to_string(colors.size()) + " colors for " +
                                    std::to_string(faces.size()) + " faces");
    // the file is float32: conversion from double only here
    std::vector<float> x(vertices.size()), y(vertices.size()), z(vertices.size());
    for (size_t i = 0; i < vertices.size(); ++i) {
        x[i] = static_cast<float>(vertices[i].x);
        y[i] = static_cast<float>(vertices[i].y);
        z[i] = static_cast<float>(vertices[i].z);
    }
    std::vector<std::vector<int>> vertex_indices(faces.size());
    std::vector<unsigned char> red(faces.size()), green(faces.size()), blue(faces.size());
    for (size_t i = 0; i < faces.size(); ++i) {
        vertex_indices[i] = {faces[i].a, faces[i].b, faces[i].c};
        red[i] = colors[i][0];
        green[i] = colors[i][1];
        blue[i] = colors[i][2];
    }
    happly::PLYData plyOut;
    plyOut.addElement("vertex", vertices.size());
    auto& v = plyOut.getElement("vertex");
    v.addProperty<float>("x", x);
    v.addProperty<float>("y", y);
    v.addProperty<float>("z", z);
    plyOut.addElement("face", faces.size());
    auto& f = plyOut.getElement("face");
    f.addListProperty<int>("vertex_indices", vertex_indices);
    f.addProperty<unsigned char>("red", red);
    f.addProperty<unsigned char>("green", green);
    f.addProperty<unsigned char>("blue", blue);
    plyOut.write(filename, happly::DataFormat::ASCII);
}
