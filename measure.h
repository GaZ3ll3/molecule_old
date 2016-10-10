/*
 * measure.h
 *
 *  Created on: Oct 9, 2016
 *      Author: Yimin Zhong
 */
#ifndef FMM_MEASURE_H
#define FMM_MEASURE_H

#include <chrono>
#include <iostream>
#include <string>
#include <iomanip>

#ifndef DISP
#define RUN(s, func){\
func;\
}
#endif


#ifdef DISP
#define RUN(s, func){ \
std::chrono::steady_clock::time_point begin =std::chrono::steady_clock::now(); \
func;\
std::chrono::steady_clock::time_point end =  end = std::chrono::steady_clock::now();\
std::cout << std::setw(15)<< s << " "  << std::setprecision(5) << std::setw(8) << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000000.0 << " seconds"<<std::endl;\
}
#endif

#endif //FMM_MEASURE_H
