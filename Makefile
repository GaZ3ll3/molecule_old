all:
	$(CXX) main.cpp -std=c++11 -fopenmp -ffast-math -O3 -I./  -lpthread -lm -ldl -o bbfmm3d
	./bbfmm3d	

clean:
	rm -f bbfmm3d
