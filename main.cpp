#include <iostream>
#include "kernel.h"

int main() {
    vector<point> s;
    vector<point> t;
    int N = 4096;
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

    double* c = new double[N * N];
    for (int i = 0; i < N * N; ++i) {
        c[i] = 1.0;
    }


    kernel bbfmm;

    bbfmm.eval = [](point& a, point& b) {
        double d = (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z)*(a.z - b.z);
        if (d < 1e-5) return 0.;
        else return 1.0/d;
        };

    bbfmm.initialize(2, s, t, &c[0], N * N , N * N, 80, 10);

    MatrixXd potentialMatrix;
    bbfmm.run(potentialMatrix);

    delete[] c;
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
