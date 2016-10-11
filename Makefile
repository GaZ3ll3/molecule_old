OMP_FLAG = -DRUN_OMP -fopenmp

all:
	$(CXX) main.cpp -std=c++11 -ffast-math -O3 -I./  -lpthread -lm -ldl -o bbfmm3d

omp:
	$(CXX) main.cpp -std=c++11 $(OMP_FLAG) -ffast-math -O3 -I./  -lpthread -lm -ldl -o bbfmm3d
clean:
	rm -f bbfmm3d
