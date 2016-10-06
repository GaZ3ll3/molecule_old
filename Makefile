
cxx = clang++-3.5
all:
	$(cxx) main.cpp tree.cpp -std=c++11 -O3 -fopenmp -I/usr/include/eigen3  -lpthread -lm -ldl -o bbfmm3d
	
