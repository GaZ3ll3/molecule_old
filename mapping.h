//
// Created by lurker on 10/24/16.
//

#ifndef BBFMM3D_STANDARD_H
#define BBFMM3D_STANDARD_H


#include <vector>
#include "node.h"
// generate standard points.

/*
 * L2 norm of a vector
 */
inline double norm(const double &x, const double &y, const double &z) noexcept {
    return sqrt(x * x + y * y + z * z);
}

inline double norm(const point &a) {
    return norm(a.x, a.y, a.z);
}

inline double inner(point a, point b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
/*
 * L2 norm squared
 */
inline double normSqr(const double &x, const double &y, const double &z) noexcept {
    return x * x + y * y + z * z;
}

inline double angle(point &a, point &b) noexcept {
    double r = norm(a.x, a.y, a.z);
    double s = norm(b.x, b.y, b.z);
    return acos((a.x * b.x + a.y * b.y + a.z * b.z) / r / s);
}

inline double theta(double psi_i, double psi_t, double phi) {
    using std::cos;
    using std::sqrt;
    using std::sin;
    double lambda = cos(psi_i) / cos(psi_t);
    double sb = sin(phi / 2.0);

    double a = sb * lambda;
    a = 4 * a * a - 1.0;
    double b = -4 * sb * sb * lambda;
    double l = (-b - sqrt(b * b - 4 * a)) / 2.0 / a;

    double _a = l * sin(psi_i);
    double _b = sqrt(4 * sb * sb - cos(psi_t) * cos(psi_t)) * l * lambda;


    double ret = (1 + _a * _a - _b * _b) / (2 * _a);
    if (abs(ret - 1.0) < 1e-12) return 0.;
    else {
        return acos(ret);
    }

}


/*
 * a single point quadrature rule for non-singular triangles
 */
inline void quadrature(vector<point> &source,
                       vector<double> &weight,
                       const double &sRadius,
                       double &triangleArea,
                       point &a,
                       point &b,
                       point &c,
                       int &triangleId) noexcept {

    double alpha, beta, gamma, psi_t, T_s;
    alpha = angle(b, c);
    beta = angle(c, a);
    gamma = angle(a, b);
    T_s = triangleArea / sRadius / sRadius;

    psi_t = asin(sqrt(cos(beta) * cos(beta) + cos(gamma) * cos(gamma) - 2 * cos(alpha) * cos(beta) * cos(gamma)) /
                 sin(alpha));

    double qqw[3] = {-sqrt(3.) / sqrt(5.), 0., sqrt(3.) / sqrt(5.)};
    double www[3] = {5. / 18., 4. / 9., 5. / 18.};
    int N = 3;
    double theta_c = 0.;
    double psi_c = 0.;

    for (int i = 0; i < N; ++i) {
        double psi_ = psi_t + (M_PI / 2.0 - psi_t) * (1.0 + qqw[i]) / 2.0;
        double theta_l = theta(psi_, psi_t, beta);
        double theta_u = alpha - theta(psi_, psi_t, gamma);
        assert(theta_u > theta_l);
        theta_c += (M_PI / 2.0 - psi_t) * (theta_u * theta_u - theta_l * theta_l) * sin(psi_) * www[i];
        psi_c += (M_PI / 2.0 - psi_t) * (theta_u - theta_l) * psi_ * sin(psi_) * www[i];
    }

    theta_c /= 2.0 * T_s;
    psi_c /= T_s;


    double mu = (sin(theta_c) * sin(psi_c) -
                 cos(psi_c) * sqrt(1 - cos(beta) * cos(beta) / sin(psi_t) / sin(psi_t)) / cos(psi_t)) / sin(alpha);
    double lambda = cos(theta_c) * sin(psi_c) - cos(beta) * cos(psi_c) / cos(psi_t) - mu * cos(alpha);
    double q = 1.0 / (lambda + mu + cos(psi_c) / cos(psi_t));
    mu *= q;
    lambda *= q;

    /*
     * barycentric coordinate of "barycentric" point
     */
//    double mu = 1.0 / 3.0;
//    double lambda = mu;

    /*
     *  Cartesian coordinate
     */
    double x_c = a.x * (1 - mu - lambda) + b.x * mu + c.x * lambda;
    double y_c = a.y * (1 - mu - lambda) + b.y * mu + c.y * lambda;
    double z_c = a.z * (1 - mu - lambda) + b.z * mu + c.z * lambda;

    /*
     * length for rescale
     */
    double r = norm(x_c, y_c, z_c);

    source.push_back(point(x_c * sRadius / r, y_c * sRadius / r, z_c * sRadius / r, triangleId));
    weight.push_back(triangleArea);
}

/*
 *  get area of spherical triangle
 *
 */
inline double area(const double &a, const double &b, const double &c, const double &sRadius) noexcept {
    double s = (a + b + c) * 0.5;
    return 4 * sRadius * sRadius *
           atan(sqrt(tan(0.5 * s) * tan(0.5 * (s - a)) * tan(0.5 * (s - b)) * tan(0.5 * (s - c))));
}


inline void zProjection(vector<point> &source,
                        vector<double> &weight,
                        vector<point> &vertices,
                        const double &fixed, const double &start, const double &delta,
                        const double &sRadius,
                        const int &N, int &curTriangleId) noexcept {

    double x, y, z, angleA, angleB, angleC, triangleArea;

    z = fixed;
    for (int i = 0; i < N; ++i) {
        x = start + i * delta;
        for (int j = 0; j < N; ++j) {
            y = start + j * delta;
            point pa(x, y, z);
            point pb(x, y + delta, z);
            point pd(x + delta, y, z);
            point pc(x + delta, y + delta, z);

            angleA = angle(pb, pc);
            angleB = angle(pc, pa);
            angleC = angle(pa, pb);
            triangleArea = area(angleA, angleB, angleC, sRadius);

            quadrature(source, weight, sRadius, triangleArea, pb, pa, pc, curTriangleId);
            curTriangleId++;
            vertices.push_back(pa);
            vertices.push_back(pb);
            vertices.push_back(pc);

            angleA = angle(pd, pc);
            angleC = angle(pd, pa);
            triangleArea = area(angleA, angleB, angleC, sRadius);

            quadrature(source, weight, sRadius, triangleArea, pd, pa, pc, curTriangleId);
            curTriangleId++;
            vertices.push_back(pa);
            vertices.push_back(pd);
            vertices.push_back(pc);
        }
    }
}

inline void xProjection(vector<point> &source,
                        vector<double> &weight,
                        vector<point> &vertices,
                        const double &fixed, const double &start, const double &delta,
                        const double &sRadius,
                        const int &N, int &curTriangleId) noexcept {

    double x, y, z, angleA, angleB, angleC, triangleArea;

    x = fixed;
    for (int i = 0; i < N; ++i) {
        y = start + i * delta;
        for (int j = 0; j < N; ++j) {
            z = start + j * delta;
            point pa(x, y, z);
            point pb(x, y + delta, z);
            point pd(x, y, z + delta);
            point pc(x, y + delta, z + delta);

            angleA = angle(pb, pc);
            angleB = angle(pc, pa);
            angleC = angle(pa, pb);

            triangleArea = area(angleA, angleB, angleC, sRadius);

            quadrature(source, weight, sRadius, triangleArea, pb, pa, pc, curTriangleId);
            curTriangleId++;
            vertices.push_back(pa);
            vertices.push_back(pb);
            vertices.push_back(pc);

            angleA = angle(pd, pc);
            angleC = angle(pd, pa);

            triangleArea = area(angleA, angleB, angleC, sRadius);

            quadrature(source, weight, sRadius, triangleArea, pd, pa, pc, curTriangleId);
            curTriangleId++;
            vertices.push_back(pa);
            vertices.push_back(pd);
            vertices.push_back(pc);
        }
    }
}

inline void yProjection(vector<point> &source,
                        vector<double> &weight,
                        vector<point> &vertices,
                        const double &fixed, const double &start, const double &delta,
                        const double &sRadius,
                        const int &N, int &curTriangleId)  noexcept {

    double x, y, z, angleA, angleB, angleC, triangleArea;

    y = fixed;
    for (int i = 0; i < N; ++i) {
        x = start + i * delta;
        for (int j = 0; j < N; ++j) {
            z = start + j * delta;
            point pa(x, y, z);
            point pb(x + delta, y, z);
            point pd(x, y, z + delta);
            point pc(x + delta, y, z + delta);

            angleA = angle(pb, pc);
            angleB = angle(pc, pa);
            angleC = angle(pa, pb);

            triangleArea = area(angleA, angleB, angleC, sRadius);

            quadrature(source, weight, sRadius, triangleArea, pb, pa, pc, curTriangleId);
            curTriangleId++;
            vertices.push_back(pa);
            vertices.push_back(pb);
            vertices.push_back(pc);

            angleA = angle(pd, pc);
            angleC = angle(pd, pa);

            triangleArea = area(angleA, angleB, angleC, sRadius);

            quadrature(source, weight, sRadius, triangleArea, pd, pa, pc, curTriangleId);
            curTriangleId++;
            vertices.push_back(pa);
            vertices.push_back(pd);
            vertices.push_back(pc);

        }
    }
}

/*
 * allocate source points and weights
 */
inline void cubeProjection(vector<point> &source,
                           vector<double> &weight,
                           vector<point> &vertices,
                           const double &cRadius, const double &sRadius,
                           const int &N) noexcept {

    /*
     * cRadius = side length of cube / 2
     */
    double fixed = cRadius;
    /*
     * size of each triangle
     */
    double delta = 2 * cRadius / N;

    /*
     * coordinate x, y, z
     */
    double x, y, z;

    /*
     * angles for spherical triangles
     */
    double angleA, angleB, angleC, triangleArea;
    /*
     * tracking triangleId
     */
    int curTriangleId = 0;
    /*
     *
     *  triangle is cut from a square.
     *
     * b |------| c
     *   |     /|
     *   |    / |
     *   |   /  |
     *   |  /   |
     *   | /    |
     *   |/     |
     * a －－－－| d
     */

    zProjection(source, weight, vertices, fixed, -fixed, delta, sRadius, N, curTriangleId);
    zProjection(source, weight, vertices, -fixed, -fixed, delta, sRadius, N, curTriangleId);
    xProjection(source, weight, vertices, fixed, -fixed, delta, sRadius, N, curTriangleId);
    xProjection(source, weight, vertices, -fixed, -fixed, delta, sRadius, N, curTriangleId);
    yProjection(source, weight, vertices, fixed, -fixed, delta, sRadius, N, curTriangleId);
    yProjection(source, weight, vertices, -fixed, -fixed, delta, sRadius, N, curTriangleId);
}

#endif //BBFMM3D_STANDARD_H
