#!/usr/bin/env python
# solve soliton BVP using spectral relaxation method

#######################################################################

import numpy as np
from numpy.linalg import solve

#######################################################################

# number of modes and compactification scale
n = 100; l = 1.4

# Chebyshev grid (for odd modes only!)
dt = (np.pi/2.0)/n; t = np.linspace(np.pi/2.0-dt/2.0, dt/2.0, n)

# construct Laplacian operator matrix
B = np.zeros([n,n])
D = np.zeros([n,n])

for i in range(0,n):
	k = 2*i+1
	B[i] = np.cos(k*t)
	D[i] = -k * (k*np.sin(t)*np.cos(k*t) + 2*np.cos(t)*np.sin(k*t)) * np.sin(t)**3/l**2

L = solve(B,D).T

#######################################################################

# field theory Lagrangian parameters
alpha = 1.0; nu = 1.0

# potential derivatives
def DV(x):
	return alpha*(x*x-nu*nu)*x

def DDV(x):
	return alpha*(3.0*x*x-nu*nu)

# initial guess
phi = B[0]

# relax the solution
for i in range(0,16):
	phi -= solve(L-np.diag(DDV(phi)), L@phi - DV(phi))

#######################################################################

import matplotlib.pyplot as plt

# sampling grid
x = l/np.tan(t)

# plot sampled solution
#plt.plot(x, phi,'o-')
plt.plot(x, phi-np.tanh(x/np.sqrt(2.0)), 'r')

plt.show()