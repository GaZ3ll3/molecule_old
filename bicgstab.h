/*
 * bicgstab.h
 *
 *  Created on: Oct 9, 2016
 *      Author: Yimin Zhong
 */

#ifndef BICGSTAB_H_
#define BICGSTAB_H_
#include <cmath>
#include <iostream>
#include "Eigen/Dense"

using namespace Eigen;
using std::sqrt;
using std::abs;

void bicgstab(std::function<VectorXd(VectorXd&)> f, VectorXd& rhs, VectorXd& x,
		const int _maxIter, const double _tol) {
    const double eps = 1e-15;
	int maxIter = _maxIter;
	double tol = _tol;

	int n = rhs.rows();
	VectorXd r = rhs - f(x);
	VectorXd r0 = r;
	double r0_sqnorm = r0.squaredNorm();
	double rhs_sqnorm = rhs.squaredNorm();

	if (abs(rhs_sqnorm) < eps) {
		x.setZero();
        return;
	}

	double rho = 1.0;
	double alpha = 1.0;
	double w = 1.0;

	VectorXd v = VectorXd::Zero(n);
	VectorXd p = VectorXd::Zero(n);
	VectorXd y(n), z(n);
	VectorXd kt(n), ks(n);
	VectorXd s(n), t(n);
	double tol2 = tol * tol;
	double eps2 = eps * eps;

	int i = 0;
	int restarts = 0;

    std::cout << "===========================" << std::endl;
    std::cout << "    iter    |  rel error   " << std::endl;

	while(r.squaredNorm()/rhs_sqnorm > tol2 && i < maxIter) {
		std::cout <<std::setw(6) << i <<  std::setw(20) << std::scientific << sqrt(r.squaredNorm()/rhs_sqnorm) << std::endl;
		double rho_old = rho;
		rho = r0.dot(r);
		if (abs(rho) < eps2 * r0_sqnorm) {
			r0 = r;
			rho = r0_sqnorm = r.squaredNorm();
			if (restarts ++ == 0) {
				i = 0;
			}
		}

		double beta = (rho/rho_old) * (alpha/w);
		p = r + beta * (p - w * v);
		y = p;
		v.noalias() = f(y);

		alpha = rho / r0.dot(v);
		s = r - alpha * v;

		z = s;
		t.noalias() = f(z);

		double tmp = t.squaredNorm();
        if (tmp > 0.) {
			w = t.dot(s)/tmp;
		}
		else {
			w = 0.;
		}

		x += alpha * y + w * z;
		r = s - w *t ;
		++i;
	}
    std::cout <<std::setw(6) << i <<  std::setw(20) << std::scientific<<sqrt(r.squaredNorm()/rhs_sqnorm) << std::endl;
    std::cout << "===========================" << std::endl;
}





#endif /* BICGSTAB_H_ */
