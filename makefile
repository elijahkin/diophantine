CC = clang++
CFLAGS = -std=c++20 -O3 -Wall -Wextra -Wpedantic -fopenmp

all: solve_omp

solve_omp: src/solve_omp.cpp
	$(CC) $(CFLAGS) -o bin/solve_omp src/solve_omp.cpp

clean:
	rm -f bin/*
