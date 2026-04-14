# Parallelization of Diophantine Solvers
In 2002, Resta and Meyrignac found the first nontrivial solutions to the Diophantine equation
```math
a_1^6 + a_2^6 = b_1^6 + b_2^6 + b_3^6 + b_4^6 + b_5^6
```
in the positive integers. Their method relied on modular restrictions for sixth powers, as well as the difference in number of terms on each side of the equation to radically shrink the search space. We notice that there is potential for parallelization in their method.
