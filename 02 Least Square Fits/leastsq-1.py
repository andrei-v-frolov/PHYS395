#!/usr/bin/env python
# generalized least square fit to supplied data (matmul version)
# run as: python leastsq.py < DATA
# expected data format: x f [sigma]

#######################################################################

# standard input access
from sys import stdin

# numerical libraries
import numpy as np
from numpy.linalg import lstsq, solve, svd, norm
from numpy.polynomial.chebyshev import chebvander, chebval

#######################################################################

# number of coefficients to fit
n = 30; epsilon = 0.0e-3

# load data from stdin
data = np.loadtxt(stdin)

# sanity check on supplied data format
pts,columns = data.shape
assert columns >= 2, f'Expecting at least 2 columns, got only {columns}'

# data to be fitted
x = data[:,0]
f = data[:,1]

# optional weights
W = 1.0/data[:,2]**2 if columns > 2 else np.ones(pts)

# expansion basis
B = chebvander(x,n-1)

#######################################################################

# NumPy has built-in solver minimizing |B*c-f|^2
#c,*r = lstsq(B,f,rcond=epsilon)

# it is better to cast the problem to nxn matrix as
A = np.matmul(B.T*W, B)
y = np.matmul(B.T*W, f)

# explicit regularization could be added, like so
#A += 1.0*np.diag(np.arange(n)**2)

# solve for best fit coefficients
c = solve(A,y)

# A could well be degenerate, use SVD to handle that
#U,S,V = svd(A, hermitian=True)
#print(f'SVD accuracy {norm(U*S@V - A)/norm(A)}')
#print(f'Condition number is {S.max()/S.min():g}')
#c = V.T*np.where(S > epsilon**2*S[0], 1.0/S, 0.0) @ U.T @ y

# this is equivalent to calling
#c,*r = lstsq(A,y,rcond=epsilon**2)

# evaluate best fit
g = chebval(x, c)

#######################################################################

import matplotlib.pyplot as plt

plt.plot(x, f)
plt.plot(x, g, linewidth=3)

# restrict x axis range
plt.xlim([-1.0,1.0])

# show in interactive console
plt.show()
