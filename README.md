# NNLS #

A non-negative least squares (NNLS) solver in Fortran. Modern examples exist for [Python](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.nnls.html) and [Julia](https://github.com/rdeits/NNLS.jl) at least (I'm sure there are others) but the original Fortran code they're based on is from the 80s and fixed format hurts me to read so I wrote this one.

compile with `make all` as usual; requires `BLAS` and `LAPACK`. Otherwise all strictly fortran 2018. Only tested on linux and with gfortran.

set `SHARED = 1` in the makefile to build as a shared library; it uses `iso_c_binding` and so should be straightforwardly callable from python using ctypes.

set `SHARED = 0` to compile standalone with `main.f` and then run `./test` to check it works.
this generates a set of 1000 matrices $A$ of random dimensions $m \times n \text{ for } m, n \in [2, 30] $ along with solution vectors $`x_{\text{ref}}`$, does $`A x_{\text{ref}} = b`$ and then runs `solve` with $A$ and $b$ to get a solution $x$.
note that there might be multiple solutions of the equations it generates; for each iteration it prints a.) $`\left|Ax - b\right|`$ and b.) $`\left|x_{\text{ref}} - x\right|`$ to stdout so you can check whether a.) the solution is reasonable and b.) whether it's the same solution or not.
any set of equations for which the maximum number of iterations is reached (currently set to $3n$ as in the scipy implementation) or the norm of $\left|Ax - b\right| > \epsilon = 1\times10^{-6}$ the program will attempt to write an output file with a load of relevant details.
