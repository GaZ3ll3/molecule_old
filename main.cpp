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
    vector<point> coarseSource;
    vector<point> coarseTarget;
    vector<point> coarseTriangle;
    vector<double> coarseWeight;

    /*
     * cube radius and sphere radius
     */
    double cRadius = 2.0;
    double sRadius = 1.5;

    /*
     * slice x slice grid on each face.
     */
    int slice = 16;

    cubeProjection(coarseSource, coarseWeight, coarseTriangle, cRadius, sRadius, slice);
    coarseTarget = coarseSource;

    vector<double> coarseNormalX(coarseSource.size());
    vector<double> coarseNormalY(coarseSource.size());
    vector<double> coarseNormalZ(coarseSource.size());

    for (int i = 0; i < coarseSource.size(); ++i) {
        coarseNormalX[i] = coarseSource[i].x / sRadius;
        coarseNormalY[i] = coarseSource[i].y / sRadius;
        coarseNormalZ[i] = coarseSource[i].z / sRadius;
    }

    int np = 3;
    int maxPoint = 160;
    int maxLevel = 10;

    double k = 0.03;
    double dE = 80.0;
    double dI = 2.0;

    /*
     * not necessary to handle singularity in function
     */
    auto eval_G0 = [&](point &a, point &b) {
        double d = (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z);
        double r = sqrt(d);
        return 1.0 / 4 / M_PI / r;
    };

    auto eval_Gk = [&](point &a, point &b) {
        double d = (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z);
        double r = sqrt(d);
        return exp(-k * r) / 4 / M_PI / r;
    };

    auto eval_pG0x = [&](point &a, point &b) {
        double d = (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z);
        double r = sqrt(d);
        return -(a.x - b.x) / 4 / M_PI / r / d;
    };

    auto eval_pG0y = [&](point &a, point &b) {
        double d = (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z);
        double r = sqrt(d);
        return -(a.y - b.y) / 4 / M_PI / r / d;
    };

    auto eval_pG0z = [&](point &a, point &b) {
        double d = (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z);
        double r = sqrt(d);
        return -(a.z - b.z) / 4 / M_PI / r / d;
    };

    auto eval_pGkx = [&](point &a, point &b) {
        double d = (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z);
        double r = sqrt(d);
        return -(a.x - b.x) * exp(-k * r) * (k * r + 1) / 4 / M_PI / r / d;
    };

    auto eval_pGky = [&](point &a, point &b) {
        double d = (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z);
        double r = sqrt(d);
        return -(a.y - b.y) * exp(-k * r) * (k * r + 1) / 4 / M_PI / r / d;
    };

    auto eval_pGkz = [&](point &a, point &b) {
        double d = (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z);
        double r = sqrt(d);
        return -(a.z - b.z) * exp(-k * r) * (k * r + 1) / 4 / M_PI / r / d;
    };


    auto Mapping = [&](vector<point> &source, vector<point> &target, vector<point> triangle, vector<double> &weight,
                       vector<double> &normalX, vector<double> &normalY, vector<double> &normalZ, VectorXd &phi) {
        int N = (int) source.size();

        assert(phi.rows() == 2 * N);

        kernel G0, Gk, pG0x, pG0y, pG0z, pGkx, pGky, pGkz;
        G0.eval = [&](point &a, point &b) {
            if (a == b) return singularIntegral(a, triangle, sRadius, eval_G0);
            else return eval_G0(a, b);
        };
        Gk.eval = [&](point &a, point &b) {
            if (a == b) return singularIntegral(a, triangle, sRadius, eval_Gk);
            else return eval_Gk(a, b);
        };
        pG0x.eval = [&](point &a, point &b) {
            if (a == b) return singularIntegral(a, triangle, sRadius, eval_pG0x);
            else return eval_pG0x(a, b);
        };
        pG0y.eval = [&](point &a, point &b) {
            if (a == b) return singularIntegral(a, triangle, sRadius, eval_pG0y);
            else return eval_pG0y(a, b);
        };
        pG0z.eval = [&](point &a, point &b) {
            if (a == b) return singularIntegral(a, triangle, sRadius, eval_pG0z);
            else return eval_pG0z(a, b);
        };
        pGkx.eval = [&](point &a, point &b) {
            if (a == b) return singularIntegral(a, triangle, sRadius, eval_pGkx);
            else return eval_pGkx(a, b);
        };
        pGky.eval = [&](point &a, point &b) {
            if (a == b) return singularIntegral(a, triangle, sRadius, eval_pGky);
            else return eval_pGky(a, b);
        };
        pGkz.eval = [&](point &a, point &b) {
            if (a == b) return singularIntegral(a, triangle, sRadius, eval_pGkz);
            else return eval_pGkz(a, b);
        };
        // f = partial phi/partial n
        // g = phi
        // Phi =  g   |  f
        VectorXd f(N);
        VectorXd gX(N), gY(N), gZ(N);
        for (int i = 0; i < N; ++i) {
            f(i) = phi(i + N) * weight[i];
            gX(i) = phi(i) * normalX[i] * weight[i];
            gY(i) = phi(i) * normalY[i] * weight[i];
            gZ(i) = phi(i) * normalZ[i] * weight[i];
        }


        G0.initialize(np, source, target, f, N, N, maxPoint, maxLevel);
        pG0x.initialize(np, source, target, gX, N, N, maxPoint, maxLevel);
        pG0y.initialize(np, source, target, gY, N, N, maxPoint, maxLevel);
        pG0z.initialize(np, source, target, gZ, N, N, maxPoint, maxLevel);

        Gk.initialize(np, source, target, f, N, N, maxPoint, maxLevel);
        pGkx.initialize(np, source, target, gX, N, N, maxPoint, maxLevel);
        pGky.initialize(np, source, target, gY, N, N, maxPoint, maxLevel);
        pGkz.initialize(np, source, target, gZ, N, N, maxPoint, maxLevel);


        VectorXd retG0(N), retGk(N), retpG0X(N), retpG0Y(N), retpG0Z(N), retpGkX(N), retpGkY(N), retpGkZ(N);
        G0.run(retG0);
        Gk.run(retGk);
        pG0x.run(retpG0X);
        pG0y.run(retpG0Y);
        pG0z.run(retpG0Z);
        pGkx.run(retpGkX);
        pGky.run(retpGkY);
        pG0z.run(retpGkZ);

        VectorXd ret(2 * N);
        ret.segment(0, N) = 0.5 * phi.segment(0, N) + retpG0X + retpG0Y + retpG0Z - retG0;
        ret.segment(N, N) = 0.5 * phi.segment(0, N) - retpGkX - retpGkY - retpGkZ + dI / dE * retGk;

        return ret;
    };


    int N = (int) coarseSource.size();
    auto coarseMap = [&](VectorXd &phi) {
        return Mapping(coarseSource, coarseTarget, coarseTriangle, coarseWeight, coarseNormalX, coarseNormalY,
                       coarseNormalZ, phi);
    };


    VectorXd input(2 * N);
    VectorXd start(2 * N);

    for (int i = 0; i < N; ++i) {
        input(i) = 1.0 / dE / 4.0 / M_PI / sRadius / (1 + k * sRadius);
        input(i + N) = -1.0 / dI / 4.0 / M_PI / sRadius / sRadius;
    }

//    start = input;
    start.setZero();

    VectorXd output(2 * N);
    output.setZero();
    for (int i = 0; i < N; ++i) {
        output(i) = 1.0 / dI / 4.0 / M_PI / sRadius;
    }

    gmres(coarseMap, output, start, 200, 40, 1e-3);

    std::cout << "coarse error :" << (start - input).norm() / input.norm() << std::endl;


}
