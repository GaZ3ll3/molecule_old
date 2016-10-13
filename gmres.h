//
// Created by lurker on 10/11/16.
//

#ifndef BBFMM3D_GMRES_H
#define BBFMM3D_GMRES_H

#include <cmath>
#include <iostream>
#include <chrono>
#include "Eigen/Dense"

using namespace Eigen;
using std::sqrt;
using std::abs;

bool gmres(std::function<VectorXd(VectorXd &)> f, VectorXd &rhs, VectorXd &x,
           const int _maxIter, const int _restart, const double _tol) {

    /*
     * for performance use
     */
    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point end;

    int maxIter = _maxIter;
    double tol = _tol;
    int iters = 0;

    int m = rhs.rows();
    VectorXd p0;
    p0 = rhs - f(x);
    VectorXd r0 = p0;

    if (abs(r0.norm() < tol)) {
        return true;
    }

    VectorXd w;
    w = VectorXd::Zero(_restart + 1);
    MatrixXd H;
    H = MatrixXd::Zero(m, _restart + 1);
    VectorXd tau;
    tau = VectorXd::Zero(_restart + 1);

    std::vector<JacobiRotation<double>> G(_restart);
    VectorXd e(m - 1);
    double beta;
    r0.makeHouseholder(e, tau.coeffRef(0), beta);
    w(0) = beta;
    H.bottomLeftCorner(m - 1, 1) = e;

    std::cout << "=======================================" << std::endl;
    std::cout << "    iter    |  abs error    |   time   " << std::endl;

    begin = std::chrono::steady_clock::now();

    for (int k = 1; k <= _restart; ++k) {
        ++iters;
        VectorXd v;
        v = VectorXd::Unit(m, k - 1);
        VectorXd workspace(m);

        for (int i = k - 1; i >= 0; --i) {
            v.tail(m - i).applyHouseholderOnTheLeft(H.col(i).tail(m - i - 1), tau.coeffRef(i), workspace.data());
        }

        VectorXd t = f(v);
        v = t;

        for (int i = 0; i < k; ++i) {
            v.tail(m - i).applyHouseholderOnTheLeft(H.col(i).tail(m - i - 1), tau.coeffRef(i), workspace.data());
        }

        if (v.tail(m - k).norm() != 0.) {
            if (k <= _restart) {
                VectorXd e(m - k - 1);
                double beta_;
                v.tail(m - k).makeHouseholder(e, tau.coeffRef(k), beta_);
                H.col(k).tail(m - k - 1) = e;
                v.tail(m - k).applyHouseholderOnTheLeft(H.col(k).tail(m - k - 1), tau.coeffRef(k), workspace.data());
            }
        }

        if (k > 1) {
            for (int i = 0; i < k - 1; ++i) {
                v.applyOnTheLeft(i, i + 1, G[i].adjoint());
            }
        }

        if (k < m && v(k) != 0.) {
            G[k - 1].makeGivens(v(k - 1), v(k));
            v.applyOnTheLeft(k - 1, k, G[k - 1].adjoint());
            w.applyOnTheLeft(k - 1, k, G[k - 1].adjoint());
        }

        H.col(k - 1).head(k) = v.head(k);

        bool stop = (k == m || abs(w(k)) < tol || iters == maxIter);
        end = std::chrono::steady_clock::now();
        std::cout << std::setw(6) << iters << std::setw(20) << std::scientific << abs(w(k))
                  << std::setw(12) << std::fixed
                  << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() / 1000000.0
                  << std::endl;
        begin = std::chrono::steady_clock::now();

        if (stop || k == _restart) {
            VectorXd y = w.head(k);
            H.topLeftCorner(k, k).template triangularView<Eigen::Upper>().solveInPlace(y);

            VectorXd x_new = y(k - 1) * VectorXd::Unit(m, k - 1);
            x_new.tail(m - k + 1).applyHouseholderOnTheLeft(H.col(k - 1).tail(m - k), tau.coeffRef(k - 1),
                                                            workspace.data());

            for (int i = k - 2; i >= 0; --i) {
                x_new += y(i) * VectorXd::Unit(m, i);
                x_new.tail(m - i).applyHouseholderOnTheLeft(H.col(i).tail(m - i - 1), tau.coeffRef(i),
                                                            workspace.data());
            }

            x += x_new;

            if (stop) {
                if (abs(w(k)) < tol) {
                    std::cout << "=============== CONVERGED =============" << std::endl;
                } else if (iters == maxIter) {
                    std::cout << "============= MAX ITERATION ==========" << std::endl;
                }
                return true;
            } else {
                k = 0;

                VectorXd p0 = f(x);
                VectorXd p1 = p0;
                r0 = rhs - p1;
                w = VectorXd::Zero(_restart + 1);
                H = MatrixXd::Zero(m, _restart + 1);
                tau = VectorXd::Zero(_restart + 1);
                double beta__;
                r0.makeHouseholder(e, tau.coeffRef(0), beta__);
                w(0) = beta__;
                H.bottomLeftCorner(m - 1, 1) = e;
            }
        }
    }
    std::cout << "====== Not Converged ======" << std::endl;
    return false;
}


#endif //BBFMM3D_GMRES_H
