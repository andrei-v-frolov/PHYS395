#!/usr/bin/env python
# compute orthogonal polynomial expansion coefficients

#######################################################################

import numpy as np
from numpy.polynomial.legendre import legvander, legval
from numpy.linalg import solve

#######################################################################

# number of coefficients
n = 30

# evaluation grid - uniform
#x = np.linspace(-1.0, 1.0, n)

# evaluation grid - cosine
theta = np.linspace(np.pi, 0.0, n)
x = np.cos(theta)

# function to approximate
f = np.exp(-x*x*4.5)

# coefficients satisfy f = B*c
B = legvander(x,n-1); c = solve(B,f)

#######################################################################

import matplotlib.pyplot as plt

# test grid
pts = 1000
x = np.linspace(-1.0, 1.0, pts)
f = np.exp(-x*x*4.5)
y = legval(x,c)

# plot residual
#plt.plot(x,f,'b-')
plt.plot(x,y-f,'r-')

# restrict x axis range
plt.xlim([-1.0,1.0])

#plt.show()

#######################################################################

from os import environ
environ['matplotlib.backend'] = 'pdf'
plt.savefig("output.pdf", bbox_inches='tight', pad_inches=0.02, transparent=True)
