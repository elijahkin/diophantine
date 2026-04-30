#!/bin/bash
#SBATCH --job-name=diophantine-omp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --exclusive
#SBATCH --time=00:10:00

set -euo pipefail

A_MAX="3000"

for t in 1 2 4 8 16 32; do
	export OMP_NUM_THREADS="$t"
	echo "Running with OMP_NUM_THREADS=$t"
	./bin/solve_omp "$A_MAX"
done
