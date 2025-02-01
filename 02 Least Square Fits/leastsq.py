#!/usr/bin/env python
# generalized least square fit to supplied data

#######################################################################

import numpy as np
from numpy.linalg import lstsq, solve, svd, norm
from numpy.polynomial.chebyshev import chebvander, chebval

#######################################################################

# number of coefficients
n = 30

# load data from file (TBD: read from stdin)
data = np.loadtxt('data.txt')

# data to be fitted (TBD: add point weights)
x = data[:,0]
f = data[:,1]

# expansion basis
B = chebvander(x,n-1)

# NumPy has built-in solver minimizing |B*c-f|^2
#c,*r = lstsq(B,f,rcond=1.0e-3)

#######################################################################

# it is better to cast the problem to nxn matrix as
A = np.matmul(B.T, B)
y = np.matmul(B.T, f)

# explicit regularization could be added, like so
#A += 1.0*np.diag(np.arange(n)**2)

# solve for best fit coefficients
c = solve(A,y)

# A could well be degenerate, use SVD to handle that
#U,S,V = svd(A, full_matrices=False, hermitian=True)
#print(f'SVD accuracy {norm(U*S@V - A)/norm(A)}')
#print(f'Condition number is {S.max()/S.min():g}')
#c = V.T*np.where(S > 1.0e-6*S[0], 1.0/S, 0.0) @ U.T @ y

# this is equivalent to calling
#c,*r = lstsq(A,y,rcond=1.0e-6)

# evaluate best fit
g = chebval(x, c)

#######################################################################

import matplotlib.pyplot as plt

plt.plot(x, f)
plt.plot(x, g)

# restrict x axis range
plt.xlim([-1.0,1.0])

# show in interactive console
plt.show()
