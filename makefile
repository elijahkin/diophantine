CC = g++
# Note: add -fopenmp for multi-threaded on Linux/Zaratan
CFLAGS = -std=c++20 -O3 -Wall -Wextra -Wpedantic

all: solve_ompORG solve_ompHARD solve_ompGEN

solve_ompORG: src/solve_ompORG.cpp
	$(CC) $(CFLAGS) -o bin/solve_ompORG src/solve_ompORG.cpp

solve_ompHARD: src/solve_ompHARD.cpp
	$(CC) $(CFLAGS) -o bin/solve_ompHARD src/solve_ompHARD.cpp

solve_ompGEN: src/solve_ompGEN.cpp
	$(CC) $(CFLAGS) -o bin/solve_ompGEN src/solve_ompGEN.cpp

clean:
	rm -f bin/*
