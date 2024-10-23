# Evaluating Wigner D-matrices

This repository contains codes in C++, Fortran and Python to 
evaluate the [Wigner D-matrix](https://en.wikipedia.org/wiki/Wigner_D-matrix).

The Fortran and Python version use the algorithm described by 
Choi et al.[^1], which resursively computes the individual matrix elements.


The C++ version uses the diagonalization approach by Feng et al.[^2]. The
relevant matrices are diagonalized using LAPACK.




[^1]: [Choi, C. H., Ivanic, J., Gordon, M. S., & Ruedenberg, K. (1999). Rapid and stable determination of rotation matrices between spherical harmonics by direct recursion. The Journal of Chemical Physics, 111(19), 8825-8831](https://doi.org/10.1063/1.480229)
[^2]: [Feng, X. M., Wang, P., Yang, W., & Jin, G. R. (2015). High-precision evaluation of Wigner's d matrix by exact diagonalization. Physical Review E, 92(4), 043307](https://link.aps.org/doi/10.1103/PhysRevE.92.043307)


