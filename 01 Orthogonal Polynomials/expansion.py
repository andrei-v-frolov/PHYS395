#!/usr/bin/env python
# compute orthogonal polynomial expansion coefficients

#######################################################################

import numpy as np
from numpy.linalg import solve, svd, norm
from numpy.polynomial.legendre import legvander, legval
from numpy.polynomial.chebyshev import chebvander, chebval

#######################################################################

# number of coefficients
n = 30

# collocation grid - uniform
#x = np.linspace(-1.0, 1.0, n)

# collocation grid - cosine
x = np.cos(np.linspace(np.pi, 0.0, n))

# function to approximate
f = np.exp(-x*x*4.5)

#######################################################################

# find coefficients that satisfy f = B*c
B = legvander(x,n-1); c = solve(B,f)

# compute singular value decomposition
#U,S,V = svd(B); print(norm(U*S@V - B)/n)

# or, just the singular values
S = svd(B, compute_uv=False)

# if S was *not* sorted already, we could do
#S = np.sort(S, kind='stable')[::-1]

# for ill-conditioned matrix, pseudo-inverse is better!
#c = V.T*np.where(S > 1.0e-12*S[0], 1.0/S, 0.0) @ U.T @ f

# compute condition number of a matrix
print(f'Condition number is {S[0]/S[-1]:.3g}')

#######################################################################

# number of test points
pts = 1024

# evaluation grid
x = np.linspace(-1.0, 1.0, pts)
f = np.exp(-x*x*4.5)

# evaluate polynomial expansion
g = legval(x,c)

#######################################################################

import matplotlib.pyplot as plt

plt.figure()

#plt.plot(x, f, '-', linewidth=3)
#plt.plot(x, g, '-', linewidth=3)
plt.plot(x, g-f, 'r-')

plt.axhline(0.0, color='black', linestyle=':', zorder=0)

# restrict x axis range
plt.xlim([-1.0,1.0])

# show in interactive console
plt.show()
