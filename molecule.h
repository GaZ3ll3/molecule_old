//
// Created by lurker on 10/24/16.
//

#ifndef BBFMM3D_MOLECULE_H
#define BBFMM3D_MOLECULE_H

#include "mapping.h"

inline void moleculeProjection(vector<point> &centers, vector<double> &radius,
                               vector<point> &source, vector<double> &weight,
                               vector<double> &NormalX, vector<double> &NormalY, vector<double> &NormalZ,
                               vector<point> &vertices, const double &cRadius, const int &N) {
    // standard unit ball with N parts.
    // for a ball with radius r, the N should change to N*r to sustain the size.
    int M = centers.size();
    int triId = 0;
    vector<vector<int>> adjacent(M);

    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < M; ++j) {
            if (j == i) continue;
            double distance = norm(centers[i].x - centers[j].x,
                                   centers[i].y - centers[j].y,
                                   centers[i].z - centers[j].z);
            if (distance < radius[i] + radius[j]) {
                adjacent[i].push_back(j);
            }
        }
    }

<<<<<<< HEAD
    std::cout << "adjacent list" << std::endl;
    for (int i = 0; i < M; ++i) {
        std::cout << i << " : ";
=======
    for (int i = 0; i < M; ++i) {
        std::cout << "atom " << i << " : ";
>>>>>>> multipole atoms exmaple
        for (int j = 0; j < adjacent[i].size(); ++j) {
            std::cout << adjacent[i][j] << " ";
        }
        std::cout << "\n";
    }


    for (int i = 0; i < M; ++i) {
        vector<point> cur_sources;
        vector<double> cur_weights;
        vector<point> cur_vertices; // each 3 forms a triangle
        double r = radius[i];
        int cur_N = int(r * N);
        // get standard points and weights at (0.,0.,0.)
        cubeProjection(cur_sources, cur_weights, cur_vertices, cRadius, r, cur_N);

        int point_N = cur_sources.size();

        for (int j = 0; j < point_N; ++j) {

            // check if this point is in or out.
            // traverse all adjacent balls.
            bool isValid = true;
            for (int k = 0; k < adjacent[i].size(); ++k) {

                if (!isValid) break;
                else {
                    int l = adjacent[i][k];
                    // vector from center to point, length is r.
                    point out(cur_sources[j].x, cur_sources[j].y, cur_sources[j].z);

                    // vector from center to center
                    point cline(centers[l].x - centers[i].x, centers[l].y - centers[i].y, centers[l].z - centers[i].z);

                    double d = norm(centers[l].x - centers[i].x, centers[l].y - centers[i].y, centers[l].z -
                                                                                              centers[i].z);

                    double theta = acos((radius[i] * radius[i] + d * d - radius[l] * radius[l]) / (2 * radius[i] * d));
                    double phi = acos(inner(out, cline) / radius[i] / d);

                    if (phi < theta) {
                        isValid = false;
                    }
                }
            }

            if (isValid) {

                NormalX.push_back((cur_sources[j].x) / radius[i]);
                NormalY.push_back((cur_sources[j].y) / radius[i]);
                NormalZ.push_back((cur_sources[j].z) / radius[i]);

                source.push_back(cur_sources[j] + centers[i]);
                source[source.size() - 1].triangleId = triId;
                triId++;
                source[source.size() - 1].ballId = i;

                weight.push_back(cur_weights[j]);

                vertices.push_back(cur_vertices[3 * j] + centers[i]);
                vertices.push_back(cur_vertices[3 * j + 1] + centers[i]);
                vertices.push_back(cur_vertices[3 * j + 2] + centers[i]);

            }
        }
    }
}


#endif //BBFMM3D_MOLECULE_H
