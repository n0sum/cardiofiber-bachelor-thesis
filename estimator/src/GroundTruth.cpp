#include "GroundTruth.h"
#include <happly.h>
#include <fstream>
#include <stdexcept>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <glm/geometric.hpp>

std::vector<FaceGroundTruth> readGroundTruthPly(const std::string& filename) {
    happly::PLYData plyIn(filename);
    auto& fe = plyIn.getElement("face");
    const std::vector<double> fx = fe.getProperty<double>("fx");
    const std::vector<double> fy = fe.getProperty<double>("fy");
    const std::vector<double> fz = fe.getProperty<double>("fz");
    const std::vector<double> cvL = fe.getProperty<double>("CV_L");
    const std::vector<double> cvT = fe.getProperty<double>("CV_T");

    std::vector<FaceGroundTruth> fGT(fx.size());
    for (int i = 0; i < fx.size(); ++i)
        fGT[i] = {glm::dvec3(fx[i], fy[i], fz[i]), cvL[i], cvT[i]};
    return fGT;
}

void exportErrorsAsCsv(const std::string& filename,
                       const cardiofiber::Mesh& mesh, const std::vector<cardiofiber::FaceAnisotropy>& est,
                       const std::vector<FaceGroundTruth>& fGT) {
    const std::vector<glm::dvec3> av = mesh.areaVectors();           // cross(e1,e2): area = |cross| / 2
    std::ofstream out(filename);
    if (!out) throw std::runtime_error("exportErrorsAsCsv: cannot open " + filename);
    out << std::scientific << std::setprecision(9);
    out << "face_id,area,estimable,cond,CV_L_est,CV_T_est,CV_L_true,CV_T_true,"
           "err_angle,err_cv_l,err_ratio\n";
    for (int i = 0; i < mesh.faces().size(); ++i) {
        const cardiofiber::FaceAnisotropy &a = est[i];
        const FaceGroundTruth &g = fGT[i];
        out << i << ',' << glm::length(av[i]) / 2.0 << ',' << (a.estimable ? 1 : 0) << ',';
        if (a.cond > 0.0) out << a.cond;                       // empty if not computed
        out << ',';
        if (a.estimable) out << a.CV_L << ',' << a.CV_T << ','; else out << ",,";
        out << g.CV_L << ',' << g.CV_T << ',';
        if (a.estimable) {                                      // empty if not estimable
            // |a.b|: angle between directions, unoriented
            const double cosAngle = std::abs(glm::dot(a.fiber_dir, g.fiber_dir));
            out << std::acos(std::min(1.0, cosAngle)) * 180.0 / M_PI << ',';
            // relative errors in %, signed, r = CV_L / CV_T; the same formulas, and the one of the angle,
            // are in estimator/roney.py: whoever changes one changes the other too
            out << 100.0 * (a.CV_L - g.CV_L) / g.CV_L << ',';
            const double rEst = a.CV_L / a.CV_T;
            const double rTrue = g.CV_L / g.CV_T;
            out << 100.0 * (rEst - rTrue) / rTrue;
        } else out << ",,";
        out << '\n';
    }
}
