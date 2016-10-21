/*
 * main.cpp
 *
 *  Created on: Oct 20, 2016
 *      Author: Yimin Zhong
 */
#ifndef BBFMM3D_EXTRAPOLATION_H
#define BBFMM3D_EXTRAPOLATION_H

#include "kernel.h"
#include "geometry.h"

inline VectorXd extrapolation(std::vector<point> &coarseSource,
                              std::vector<point> &fineSource,
                              std::vector<point> &target,
                              std::vector<point> &coarseTriangle,
                              std::vector<point> &fineTriangle,
                              std::vector<double> &coarseWeight,
                              std::vector<double> &fineWeight,
                              std::function<double(point &, point &)> eval,
                              const double &sRadius, int np, int maxPoint, int maxLevel) noexcept {

    kernel G_coarse, G_fine;
    G_coarse.eval = [&](point &a, point &b) {
        if (a == b) // a = b
        {
            return singularIntegral(a, coarseTriangle, sRadius, eval);
        } else return eval(a, b);
    };

    G_fine.eval = [&](point &a, point &b) {
        if (a == b) // a = b
        {
            return singularIntegral(a, fineTriangle, sRadius, eval);
        } else return eval(a, b);
    };


    int N_coarse = (int) coarseWeight.size();
    int N_fine = (int) fineWeight.size();
    VectorXd f(coarseWeight.size());
    VectorXd g(fineWeight.size());

    for (int i = 0; i < N_coarse; ++i) {
        f(i) = coarseWeight[i];
    }

    for (int i = 0; i < N_fine; ++i) {
        g(i) = fineWeight[i];
    }

    G_coarse.initialize(np, coarseSource, target, f, N_coarse, N_coarse, maxPoint, maxLevel);
    G_fine.initialize(np, fineSource, target, g, N_fine, N_coarse, maxPoint, maxLevel);
    VectorXd potential;
    VectorXd potential2;
    G_coarse.run(potential);
    G_fine.run(potential2);
    return potential2 * 2.0 - potential;

}

#endif //BBFMM3D_EXTRAPOLATION_H
