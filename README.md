# optimization-project

### Optimization Algorithm Benchmarking: A Case Study on the Bohachevsky 1 Problem (BF1)

### Purpose of repository
This repository created for benchmarking some optimization algorithm on Bohachevsky 1 Problem (BF1).

**Available optimization algorithm:**

- Newton-Raphson
- Hestenes-Stiefel
- Polak-Ribiere
- Fletcher-Reeves
- Gradient-Descent

**Helper Functions:**

- func.m
- gradfunc.m
- hessianfunc.m

### Usage

The Matlab script files given below can be run directly in Matlab or any other proper environment.

**Matlab script files:**

- newton_raphson.m
- hestenes_stiefel.m
- polak_ribiere.m
- fletcher_reeves.m
- gradient_descent.m

**Output:**

In these script files,
- the success rate of the algorithm (based on the maximum number of iterations, which can be configured within the script),
- the average number of iterations for each initial value,
- the average elapsed time for each initial value

are calculated. Upon completion of the script execution, these results are printed to the command window.
