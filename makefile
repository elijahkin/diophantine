CC = clang++
CFLAGS = -std=c++20 -O3 -Wall -Wextra -Wpedantic -fopenmp

all: diophantine

diophantine: src/diophantine.cpp
	$(CC) $(CFLAGS) -o bin/diophantine src/diophantine.cpp

clean:
	rm -f bin/*
