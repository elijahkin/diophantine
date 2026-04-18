# Parallelization of Diophantine Solvers
In 2002, Resta and Meyrignac found the first nontrivial solutions to the Diophantine equation
```math
a_1^6 + a_2^6 = b_1^6 + b_2^6 + b_3^6 + b_4^6 + b_5^6
```
in the positive integers. Their method relied on modular restrictions for sixth powers, as well as the difference in number of terms on each side of the equation to radically shrink the search space. We notice that there is potential for parallelization in their method, as we have independent iterations for each possible value of $a_1$.

For symmetry breaking, our solver enforces $a_1 \geq a_2$, as well as $b_1 \geq b_2$ and $b_3 \geq b_4 \geq b_5$.
## Usage
The solver accepts one command line argument indicating the maximum value of $a_1$. For example, if we run `./bin/solve_omp 3000`, we obtain the following result:
```
1117^6 + 770^6 = 602^6 + 212^6 + 1092^6 + 861^6 + 84^6
1117^6 + 770^6 = 861^6 + 212^6 + 1092^6 + 602^6 + 84^6
1117^6 + 770^6 = 1092^6 + 212^6 + 861^6 + 602^6 + 84^6
2041^6 + 691^6 = 1893^6 + 1468^6 + 1407^6 + 1302^6 + 1246^6
2441^6 + 752^6 = 2096^6 + 1266^6 + 2184^6 + 1484^6 + 1239^6
2827^6 + 151^6 = 1488^6 + 390^6 + 2653^6 + 2296^6 + 1281^6
2959^6 + 2470^6 = 2481^6 + 850^6 + 2954^6 + 798^6 + 420^6
Finished search in 10824 ms!
```
