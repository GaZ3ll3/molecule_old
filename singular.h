//
// Created by lurker on 10/24/16.
//

#ifndef BBFMM3D_SINGULAR_H
#define BBFMM3D_SINGULAR_H

#include "mapping.h"

inline double singularAngle(point &a, point &b, int id, vector<point> &centers) {
    point pa = a - centers[id];
    point pb = b - centers[id];
    double r = norm(pa);
    double s = norm(pb);
    return acos(inner(pa, pb) / r / s);
}


inline void gauss(point &singularPoint, point &left, point &right,
                  const double &sRadius, vector<point> &centers,
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

            int id = singularPoint.ballId;
            point pp(px, py, pz);
            point pq = pp - centers[id];

            /*
             * project to sphere
             */
            double scale = sRadius / norm(pq);
            quadraturePoint.push_back(point(pq.x * scale, pq.y * scale, pq.z * scale) + centers[id]);
        }
    }
}


/*
 * eval should be in form of eval(source, target)
 */
inline double singularIntegral(point &singular, vector<point> &vertices, const double &sRadius, vector<point> &centers,
                               std::function<double(point &, point &)> eval) noexcept {
    /*
     * get 3 vertices of triangle enclosing singularity
     */
    point &pointA = vertices[3 * singular.triangleId];
    point &pointB = vertices[3 * singular.triangleId + 1];
    point &pointC = vertices[3 * singular.triangleId + 2];

    int id = singular.ballId;

    double angleA = singularAngle(pointB, pointC, id, centers);
    double angleB = singularAngle(pointC, pointA, id, centers);
    double angleC = singularAngle(pointA, pointB, id, centers);
    double triangleArea = area(angleA, angleB, angleC, sRadius);

    double angleSubA, angleSubB, angleSubC, triangleSubArea;

    angleSubA = singularAngle(pointA, singular, id, centers);
    angleSubB = singularAngle(pointB, singular, id, centers);
    angleSubC = singularAngle(pointC, singular, id, centers);

    vector<point> quadraturePoint;
    vector<double> quadratureWeight;

    double ret = 0.;

    /*
     * AB Part
     */
    gauss(singular, pointA, pointB, sRadius, centers, quadraturePoint, quadratureWeight);
    triangleSubArea = area(angleC, angleSubA, angleSubB, sRadius);

    for (int i = 0; i < quadraturePoint.size(); ++i) {
        ret += quadratureWeight[i] * triangleSubArea / triangleArea * eval(quadraturePoint[i], singular);
    }
    quadraturePoint.clear();
    quadratureWeight.clear();

    /*
     * BC part
     */
    gauss(singular, pointB, pointC, sRadius, centers, quadraturePoint, quadratureWeight);
    triangleSubArea = area(angleA, angleSubB, angleSubC, sRadius);

    for (int i = 0; i < quadraturePoint.size(); ++i) {
        ret += quadratureWeight[i] * triangleSubArea / triangleArea * eval(quadraturePoint[i], singular);
    }
    quadraturePoint.clear();
    quadratureWeight.clear();


    /*
     * CA part
     */
    gauss(singular, pointC, pointA, sRadius, centers, quadraturePoint, quadratureWeight);
    triangleSubArea = area(angleB, angleSubC, angleSubA, sRadius);

    for (int i = 0; i < quadraturePoint.size(); ++i) {
        ret += quadratureWeight[i] * triangleSubArea / triangleArea * eval(quadraturePoint[i], singular);
    }
    quadraturePoint.clear();
    quadratureWeight.clear();

    return ret;
}


#endif //BBFMM3D_SINGULAR_H
