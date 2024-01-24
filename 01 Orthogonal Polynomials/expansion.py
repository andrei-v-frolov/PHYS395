#!/usr/bin/env python
# compute orthogonal polynomial expansion coefficients

import numpy as np
from numpy.polynomial.legendre import legvander, legval
from numpy.linalg import solve

# number of coefficients
n = 30

# evaluation grid - uniform
x = np.linspace(-1.0, 1.0, n)

# evaluation grid - Chebyshev
theta = np.linspace(np.pi, 0.0, n)
x = np.cos(theta)

# function to approximate
f = np.exp(-x*x*4.5)

# coefficients satisfy A*c = f
A = legvander(x,n-1)
c = solve(A,f)

#######################################################################

import matplotlib.pyplot as plt

#plt.plot(x,f,'g.')

# test grid
pts = 1000
x = np.linspace(-1.0, 1.0, pts)
f = np.exp(-x*x*4.5)
y = legval(x,c)

# plot residual
#plt.plot(x,f,'b-')
plt.plot(x,y-f,'r-')

plt.show()
