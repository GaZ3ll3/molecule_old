/*
 * main.cpp
 *
 *  Created on: Oct 9, 2016
 *      Author: Yimin Zhong
 */
#include <iostream>
#include "kernel.h"
#include "bicgstab.h"
#include "gmres.h"
#include "geometry.h"
#include "extrapolation.h"

int main() {
#ifdef RUN_OMP
    omp_set_num_threads(omp_get_max_threads());
#endif
    /*
     * allocating source, target locations and weights
     */
    vector<point> coarseSource, fineSource, target, coarseTriangle, fineTriangle;
    vector<double> coarseWeight, fineWeight;

    /*
     * cube radius and sphere radius
     */
    double cRadius = 1.0;
    double sRadius = 1.0;

    /*
     * slice x slice grid on each face.
     */
    int slice = 32;

    cubeProjection(coarseSource, coarseWeight, coarseTriangle, cRadius, sRadius, slice);
    cubeProjection(fineSource, fineWeight, fineTriangle, cRadius, sRadius, 2 * slice);

    target = coarseSource;

    int np = 4;
    int maxPoint = 160;
    int maxLevel = 10;

    /*
     * not necessary to handle singularity in function
     */
    auto eval = [&](point &a, point &b) {
        double d = (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z);
        double r = sqrt(d);
        return 1.0 / 4 / M_PI / r;
    };

    VectorXd ret = extrapolation(coarseSource, fineSource, target, coarseTriangle, fineTriangle, coarseWeight,
                                 fineWeight, eval, sRadius, np, maxPoint, maxLevel);

    std::cout << (ret - VectorXd::Ones(coarseSource.size())).norm() / std::sqrt(coarseSource.size()) << std::endl;
}
