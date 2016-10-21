/*
 * geometry.h
 *
 *  Created on: Oct 9, 2016
 *  Updated on: Oct 20, 2016
 *      Author: Yimin Zhong
 */

#ifndef FMM_GEOMETRY_H
#define FMM_GEOMETRY_H

#include <vector>
#include "node.h"

/*
 * L2 norm of a vector
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
 * angle between two vectors
 */
inline double angle(point &a, point &b) noexcept {
    double r = norm(a.x, a.y, a.z);
    double s = norm(b.x, b.y, b.z);
    return acos((a.x * b.x + a.y * b.y + a.z * b.z) / r / s);
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

    /*
     * barycentric coordinate of "barycentric" point
     */
    double mu = 1.0 / 3.0;
    double lambda = mu;

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

inline void gauss(point &singularPoint, point &left, point &right,
                  const double &sRadius,
                  vector<point> &quadraturePoint, vector<double> &quadratureWeight) noexcept {
    /*
     * high order accuracy scheme with 3 point each dimension.
     */
    vector<double> gaussPoint = {(1 - std::sqrt(3.) / sqrt(5.)) / 2., 0.5, (1 + std::sqrt(3.) / sqrt(5.)) / 2.};
    vector<double> gaussWeight = {5. / 18., 4. / 9., 5. / 18.};

    assert(gaussWeight.size() == gaussPoint.size());

    int N = (int) gaussPoint.size();

    // Cartesian tensor of gauss points
    for (int i = 0; i < N; ++i) {
        double curWeight = gaussWeight[i];
        double curPoint1 = gaussPoint[i];
        for (int j = 0; j < N; ++j) {
            double combinedWeight = gaussWeight[j] * curWeight;
            /*
             * rescale due to singularity
             */
            double curPoint2 = gaussPoint[j] * (1.0 - curPoint1);
            quadratureWeight.push_back(2.0 * (1.0 - curPoint1) * combinedWeight);

            /*
             * Cartesian coordinate of quadrature point
             */
            double px = singularPoint.x * curPoint1 + right.x * curPoint2 + left.x * (1.0 - curPoint1 - curPoint2);
            double py = singularPoint.y * curPoint1 + right.y * curPoint2 + left.y * (1.0 - curPoint1 - curPoint2);
            double pz = singularPoint.z * curPoint1 + right.z * curPoint2 + left.z * (1.0 - curPoint1 - curPoint2);
            /*
             * project to sphere
             */
            double scale = sRadius / norm(px, py, pz);
            quadraturePoint.push_back(point(px * scale, py * scale, pz * scale));
        }
    }
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

/*
 * eval should be in form of eval(source, target)
 */
inline double singularIntegral(point &singular, vector<point> &vertices, const double &sRadius,
                               std::function<double(point &, point &)> eval) {
    /*
     * get 3 vertices of triangle enclosing singularity
     */
    point &pointA = vertices[3 * singular.triangleId];
    point &pointB = vertices[3 * singular.triangleId + 1];
    point &pointC = vertices[3 * singular.triangleId + 2];

    double angleA = angle(pointB, pointC);
    double angleB = angle(pointC, pointA);
    double angleC = angle(pointA, pointB);
    double triangleArea = area(angleA, angleB, angleC, sRadius);

    double angleSubA, angleSubB, angleSubC, triangleSubArea;

    angleSubA = angle(pointA, singular);
    angleSubB = angle(pointB, singular);
    angleSubC = angle(pointC, singular);

    vector<point> quadraturePoint;
    vector<double> quadratureWeight;

    double ret = 0.;

    /*
     * AB Part
     */
    gauss(singular, pointA, pointB, sRadius, quadraturePoint, quadratureWeight);
    triangleSubArea = area(angleC, angleSubA, angleSubB, sRadius);

    for (int i = 0; i < quadraturePoint.size(); ++i) {
        ret += quadratureWeight[i] * triangleSubArea / triangleArea * eval(quadraturePoint[i], singular);
    }
    quadraturePoint.clear();
    quadratureWeight.clear();

    /*
     * BC part
     */
    gauss(singular, pointB, pointC, sRadius, quadraturePoint, quadratureWeight);
    triangleSubArea = area(angleA, angleSubB, angleSubC, sRadius);

    for (int i = 0; i < quadraturePoint.size(); ++i) {
        ret += quadratureWeight[i] * triangleSubArea / triangleArea * eval(quadraturePoint[i], singular);
    }
    quadraturePoint.clear();
    quadratureWeight.clear();


    /*
     * CA part
     */
    gauss(singular, pointC, pointA, sRadius, quadraturePoint, quadratureWeight);
    triangleSubArea = area(angleB, angleSubC, angleSubA, sRadius);

    for (int i = 0; i < quadraturePoint.size(); ++i) {
        ret += quadratureWeight[i] * triangleSubArea / triangleArea * eval(quadraturePoint[i], singular);
    }
    quadraturePoint.clear();
    quadratureWeight.clear();

    return ret;
}


#endif //FMM_GEOMETRY_H
