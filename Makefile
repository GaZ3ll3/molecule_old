all:
	$(CXX) main.cpp -std=c++11 -O3 -I./  -lpthread -lm -ldl -o bbfmm3d
	./bbfmm3d	
