
cxx = clang++
all:
	$(cxx) main.cpp tree.cpp -std=c++11 -O3 -I./  -lpthread -lm -ldl -o bbfmm3d
	
