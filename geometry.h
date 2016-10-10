/*
 * geometry.h
 *
 *  Created on: Oct 9, 2016
 *      Author: Yimin Zhong
 */

#ifndef FMM_GEOMETRY_H
#define FMM_GEOMETRY_H

#include <vector>
#include "node.h"

/*
 * weight function
 */
inline double wcos(const double r) {
    if (abs(r) >= 1.0) return 0.;
    else {
        return (1.0 + cos(M_PI * r)) / 2.0;
    }
}

/*
 * L2 norm of vector
 */
inline double norm(const double &x, const double &y, const double &z) noexcept {
    return sqrt(x * x + y * y + z * z);
}

/*
 * L2 norm squared
 */
inline double normSqr(const double &x, const double &y, const double &z) noexcept {
    return x * x + y * y + z * z;
}

/*
 * calculate Jacobian, eta is norm of a point(vector).
 */
inline double Jacobian(const double radius, const double eta) {
    double kappa = 1.0 - eta / radius;
    return kappa * kappa;
}

/*
 * select points in tube near sphere.
 * todo: slow selection currently.
 */
void
spherical_tube(vector<point> &source, vector<double> &weight, double thickness, double radius, double dx) noexcept {
    double inner_radius = radius - thickness;
    double outer_radius = radius + thickness;

    double dx_cube = dx * dx * dx / thickness;

    int neg_bound = -outer_radius / dx;
    int pos_bound = -neg_bound + 1;

    double x, y, z;
    for (int i = neg_bound; i < pos_bound; ++i) {
        x = i * dx;
        for (int j = neg_bound; j < pos_bound; ++j) {
            y = j * dx;
            for (int k = neg_bound; k < pos_bound; ++k) {
                z = k * dx;
                double r2 = normSqr(x, y, z);
                if (r2 < outer_radius * outer_radius && r2 > inner_radius * inner_radius) {
                    // projection
                    double norm_p = norm(x, y, z);
                    assert(norm_p > __eps);
                    double cur_weight =
                            Jacobian(radius, norm_p - radius) * wcos((norm_p - radius) / thickness) * dx_cube;
                    // setup weight
                    weight.push_back(cur_weight);
                    // setup source
                    source.push_back(point(
                            x * radius / norm_p,
                            y * radius / norm_p,
                            z * radius / norm_p));
                }
            }
        }
    }
}


#endif //FMM_GEOMETRY_H
