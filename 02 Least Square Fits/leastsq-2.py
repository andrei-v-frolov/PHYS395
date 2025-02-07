#!/usr/bin/env python
# generalized least square fit to supplied data (serial version)
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

# accumulators
A = np.zeros([n,n])
y = np.zeros(n)

# process data line by line (NOT loading to memory)
for line in stdin:
	# parse whitespace separated floats
	x,f,*sigma = [float(s) for s in line.split()]
	# assign weight if sigma was supplied
	w = 1.0/sigma[0]**2 if len(sigma) > 0 else 1.0
	# basis is Chebyshev polynomials
	b = np.cos(np.arange(n)*np.arccos(x))
	# accumulate fit matrices
	A += w*np.outer(b,b); y += w*f*b

# explicit regularization could be added, like so
#A += 1.0*np.diag(np.arange(n)**2)

#######################################################################

# solve for best fit coefficients
c = solve(A,y)

# A could well be degenerate, use SVD to handle that
#U,S,V = svd(A, hermitian=True)
#print(f'SVD accuracy {norm(U*S@V - A)/norm(A)}')
#print(f'Condition number is {S.max()/S.min():g}')
#c = V.T*np.where(S > epsilon**2*S[0], 1.0/S, 0.0) @ U.T @ y

# this is equivalent to calling
#c,*r = lstsq(A,y,rcond=epsilon**2)

# test grid (we did NOT store original data)
x = np.linspace(-1.0,1.0,256)

# evaluate best fit 
g = chebval(x, c)

#######################################################################

import matplotlib.pyplot as plt

plt.plot(x, g, linewidth=3)

# restrict x axis range
plt.xlim([-1.0,1.0])

# show in interactive console
plt.show()
