cmake_minimum_required(VERSION 2.8)
project(bbfmm3d)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_RELEASE} -std=c++11")

set(SOURCE_FILES main.cpp geometry.h)
add_executable(bbfmm3d ${SOURCE_FILES})