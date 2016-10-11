/*
 * main.cpp
 *
 *  Created on: Oct 9, 2016
 *      Author: Yimin Zhong
 */
#include <iostream>
#include "kernel.h"
#include "bicgstab.h"
#include "geometry.h"
int main() {
#ifdef RUN_OMP
    omp_set_num_threads(omp_get_max_threads());
#endif
    /*
     * allocating source, target locations and weights
     */
    vector<point> source, target;
    vector<double> weight;

    /*
     * select points inside tube and assign to source and target.
     *
     * 1. source points = target points for simplification.
     * 2. weight is calculated by using wcos.
     * 3. thickness of tube is 0.1
     * 4. radius of sphere is 1.0
     * 5. dx=dy=dz are set to be 1.0/40.0
     */
    double thickness = 0.1;
    double radius = 1.0;
    double dx = 1.0 / 40.0;

    spherical_tube(source, weight, thickness, radius, dx);
    target = source;

    /*
     * number of points(N), and number Cheyshev points (np) used for approximation.
     * tau is the parameter for K_{\tau}.
     */
    int N = (int) source.size();
    int np = 3;
    double tau = dx;

    /*
     * allocate rhs vector for solving following second kind integral equation
     *
     *  Ax = (0.5 I + K) x = c
     *
     * here c is set to be all ones. In theory, x should be all ones.
     */
    VectorXd c(N);
    c.setOnes();


    /*
     * setting up A operator as f in above context.
     *
     * where K is set to be approximated with FMM. For detail, see below.
     */
    auto f = [&](VectorXd& v) {
        kernel bbfmm;

        /*
         * Eval function is
         *
         * K(x, y) = - r.n/r^3/(4\pi)
         *
         *  where r = |x-y|. n = n_y.
         */
        bbfmm.eval = [&](point &a, point &b) {
            double d = (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z)*(a.z - b.z);
            double r = sqrt(d);
            if (r < tau) return (1.0 / 8.0 / M_PI / tau) - (5.0 / 256.0 + 25.0 / 768.0) * tau / M_PI;
            else return -((a.x - b.x) * b.x + (a.y - b.y) * b.y + (a.z - b.z) * b.z) / d / r / (4 * M_PI);
            };

        VectorXd weighted_f(N);
        for (int i = 0; i < N; ++i) {
            weighted_f(i) = v(i) * weight[i];
        }

        bbfmm.initialize(np, source, target, weighted_f, N, N, 80, 10);

        VectorXd potentialMatrix;

        bbfmm.run(potentialMatrix);

        /*
         * plus 0.5 Identity.
         */
        potentialMatrix += v / 2.0;

        return potentialMatrix;
    };

    /*
     * allocation for solution x
     */
    VectorXd x(N);

    /*
     * using BICGSTAB to solve x.
     */
    bicgstab(f, c, x, 100, 1e-12);

    /*
     * check if norm(x - 1.0) is small enough.
     */
    VectorXd true_sol(N);
    true_sol.setOnes();
    std::cout << "BICGSTAB solution's relative error is " << std::scientific << (x - true_sol).norm() / true_sol.norm()
              << std::endl;


    /*
     *  Check solution against random target points.
     *
     *  \int K(x, y) dy = 0.  for inside
     *                  = 0.5 on edge
     *                  = 1.0 for outside
     *
     *
     */
    kernel check_bbfmm;
    check_bbfmm.eval = [&](point &a, point &b) {
        double d = (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z);
        double r = sqrt(d);
        if (r < tau) return (1.0 / 8.0 / M_PI / tau) - (5.0 / 256.0 + 25.0 / 768.0) * tau / M_PI;
        else return -((a.x - b.x) * b.x + (a.y - b.y) * b.y + (a.z - b.z) * b.z) / d / r / (4 * M_PI);
    };

    /*
     * setting up check_target for testing.
     *
     * multiply a scalar(0.9) in order to prevent points out of box [-1, 1]^3.
     */
    vector<point> check_target;
    int M = 200;
    for (int i = 0; i < M; ++i) {
        check_target.push_back(point(
                1.8 * rand() / RAND_MAX - 0.9,
                1.8 * rand() / RAND_MAX - 0.9,
                1.8 * rand() / RAND_MAX - 0.9));
    }

    VectorXd wf(N);
    for (int i = 0; i < N; ++i) {
        wf(i) = x(i) * weight[i];
    }

    /*
     * calculate \int K(x, y) dy for all x \in check_target.
     */
    check_bbfmm.initialize(np, source, check_target, wf, N, M, 80, 10);
    VectorXd ret;
    check_bbfmm.run(ret);


    /*
     * direct summation for checking
     */
    VectorXd check_ret(check_target.size());
    check_ret.setZero();
    for (int j = 0; j < check_ret.rows(); ++j) {
        for (int i = 0; i < N; ++i) {
            check_ret(j) += check_bbfmm.eval(source[i], check_target[j]) * x(i) * weight[i];
        }
    }

    std::cout << "L∞ error of above test is " << std::scientific << (ret - check_ret).lpNorm<Infinity>() << std::endl;
}
