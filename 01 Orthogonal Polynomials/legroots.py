#!/usr/bin/env python
# compute roots of Legendre polynomial using Jacobi matrix

#######################################################################

import numpy as np
from scipy.linalg import eigh_tridiagonal, norm
from numpy.polynomial.legendre import legval

#######################################################################

# order of a polynomial
n = 30

# monic Legendre polynomial satisfies recursion relation
# P[n+1] = x P[n] - n*n/(4*n*n-1) P[n-1]
k = np.arange(1,n); q = k/np.sqrt(4.0*k*k-1.0)

# find eigenvalues of (tridiagonal) Jacobi matrix
x = eigh_tridiagonal(np.zeros(n), q, eigvals_only=True)

# evaluate Legendre polynomial on found roots
c = np.zeros(n+1); c[-1] = 1.0; P = legval(x,c)

print(f'Norm of P[{n}] evaluated on found roots: {norm(P)/n}')
