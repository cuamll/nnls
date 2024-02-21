# NNLS #

A non-negative least squares (NNLS) solver in fortran.

compile with `make all` as usual; requires `BLAS` and `LAPACK`, otherwise all strictly fortran 2018.

set `SHARED = 1` in the makefile to build as a shared library; it uses `iso_c_binding` and so should be straightforwardly callable from python using ctypes.

set `SHARED = 0` to compile standalone with `main.f` and then run `./test` to check it works.
this will generate a set of 1000 matrices $A$ and solution vectors $`x_{\text{ref}}`$, does $`A x_{\text{ref}} = b`$ and then runs `solve` with $A$ and $b$ to get a solution $x$.
note that there might be multiple solutions of the equations it generates; for each iteration it prints $`\left|Ax -b\right|`$ and $`\left|x_{\text{ref}} - x\right|`$ so you can check the solution is reasonable and whether it's the same solution or not.
