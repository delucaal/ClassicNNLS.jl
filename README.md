# ClassicNNLS.jl: Non-Negative Least Squares for Julia

This package implements the "classic" non-negative least squares solver from Lawson and Hanson [1]. 
It is currently developed and tested for Julia 1.0.

Usage: 

using ClassicNNLS

ClassicNNLSSolver(X,Y,MAX_ITER=1000,tol=1e-8)

where X is the problem matrix, and Y the observed values. MAX_ITER is the maximum number of iterations
allowed in the outer cycle, tol is the maximum allowed difference between observations and predictions.

# References

[1] Lawson, C.L. and R.J. Hanson, Solving Least-Squares Problems, Prentice-Hall, Chapter 23, p. 161, 1974