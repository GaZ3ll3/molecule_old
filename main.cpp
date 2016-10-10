#include <iostream>

#include "kernel.h"
#include "bicgstab.h"

int main() {
    vector<point> s;
    vector<point> t;
    int N = 64;
    for (int i = 0; i < N; ++i) {
        double phi = M_PI * (double)rand() / RAND_MAX;
        for (int j = 0; j < N; ++j) {
            double theta = 2 * M_PI * (double)rand() / RAND_MAX;
            s.push_back(point(cos(theta) * sin(phi),
                              sin(theta) * sin(phi),
                              cos(phi)));
        }
    }
    for (int i = 0; i < N; ++i) {
        double phi = M_PI * (double)rand() / RAND_MAX;
        for (int j = 0; j < N; ++j) {
            double theta = 2 * M_PI * (double)rand() / RAND_MAX;
            t.push_back(point(cos(theta) * sin(phi),
                              sin(theta) * sin(phi),
                              cos(phi)));
        }
    }

    VectorXd c = VectorXd::Ones(N * N);




//    bbfmm.initialize(2, s, t, c, N * N , N * N, 80, 10);
//    VectorXd potentialMatrix;
//    bbfmm.run(potentialMatrix);


    auto f = [&](VectorXd& v) {
        kernel bbfmm;

        bbfmm.eval = [](point& a, point& b) {
            double d = (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z)*(a.z - b.z);
            double r = sqrt(d);
            if (d < 1e-12) return 0.;
            else return -((a.x - b.x) * b.x + (a.y - b.y)*b.y + (a.z - b.z) * b.z)/d/r;
            };

        bbfmm.initialize(4, s, t, v, N * N , N * N, 80, 10);
        VectorXd potentialMatrix;
        bbfmm.run(potentialMatrix);
        potentialMatrix /= (N*N);
        potentialMatrix +=  0.5 * v;
        return potentialMatrix;
    };

    VectorXd x(N * N);
    bicgstab(f, c, x, 100, 1e-10);
//    std::cout << x << std::endl;

//
//    kernel bbfmm2;
//    bbfmm2.eval = [](point& a, point& b) {
//        double d = (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z)*(a.z - b.z);
//        double r = sqrt(d);
//        if (d < 1e-12) return 0.;
//        else return -((a.x - b.x) * b.x + (a.y - b.y)*b.y + (a.z - b.z) * b.z)/d/r;
//    };
//
//    vector<point> tt;
//    tt.push_back(point(0., 0., 0.));
//    tt.push_back(point(-1., -1., 1.));
//    tt.push_back(point(0.3, 0.4, 0.));
//    bbfmm2.initialize(4, s, tt, x, N * N , 3, 80, 10);
//    VectorXd ret;
//    bbfmm2.run(ret);
//    std::cout << ret/(N*N) << std::endl;
//
//    double ss = 0.;
//    for (int i = 0; i < N*N;++i) {
//        ss += bbfmm2.eval(s[i], tt[1]) * x(i);
//    }
//    std::cout << ss/(N*N) << std::endl;





/*
    MatrixXd potentialMatrix2;
    potentialMatrix2 = MatrixXd::Zero(N*N, 1);

    for (int i = 0; i < N*N; ++i) {
        for (int j = 0; j < N*N; ++j) {
            potentialMatrix2(j, 0) += bbfmm.eval(s[i], t[j]);
        }
    }

    std::cout << (potentialMatrix - potentialMatrix2).norm()/potentialMatrix.norm() << std::endl;
*/
}
