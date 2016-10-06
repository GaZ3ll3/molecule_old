//
// Created by lurker on 10/2/16.
//

#ifndef FMM_MEASURE_H
#define FMM_MEASURE_H

#include <chrono>
#include <iostream>
#include <string>
#include <iomanip>

#define RUN(s, func) tic(s); func;toc();

using namespace std;

class measure {
private:
    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point end;
    std::string title;
public:
    void tic(std::string s) {
        begin = std::chrono::steady_clock::now();
        title = s;
    }

    void toc() {
        end = std::chrono::steady_clock::now();
        std::cout << setw(15)<< title << " "  << setprecision(5) << setw(8) << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000000.0
                  << " seconds"<<std::endl;
    }
};

#endif //FMM_MEASURE_H
