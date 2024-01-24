#!/usr/bin/env python
# compute orthogonal polynomials by recursion

import numpy as np

# order to expand to and number of points
nmax = 10; pts = 101

# evaluation grid
x = np.linspace(-1.0, 1.0, pts)

# initialize the array
P = np.zeros([nmax,pts])

# initialize the recursion relation
P[0,:] = 1.0
P[1,:] = x

# Chebyshev: T[n+1] = 2x*T[n] - T[n-1]
for n in range(1,nmax-1):
	P[n+1,:] = 2.0*x*P[n,:] - P[n-1,:]

# print the result
print(P)


#######################################################################

import matplotlib.pyplot as plt

for n in range(0,nmax):
	plt.plot(x,P[n,:])

plt.show()
