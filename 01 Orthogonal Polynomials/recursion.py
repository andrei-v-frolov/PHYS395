#!/usr/bin/env python
# compute orthogonal polynomials by recursion

#######################################################################

import numpy as np
from time import perf_counter as now

# order to expand to and number of points
n = 5; pts = 1024

# evaluation grid
x = np.linspace(-1.0, 1.0, pts)

# allocate the array
#P = np.empty([n+1,pts])

# it is safer to initialize it!
P = np.zeros([n+1,pts])

# time the computation
t1 = now()

# initialize the recursion relation
P[0] = 1.0
P[1] = x

# recursion for Chebyshev and Legendre polynomials
for k in range(1,n):
	#P[k+1] = 2.0*x*P[k] - P[k-1]
	P[k+1] = (2*k+1)/(k+1)*x*P[k] - k/(k+1)*P[k-1]

# wall time elapsed (note formatted string literal!)
t2 = now()
print(f'Recursion takes {t2-t1}s to complete')

#######################################################################

from numpy.polynomial.legendre import legvander

# time the computation
t1 = now()

# direct evaluation of Chebyshev polynomials
#theta = np.arccos(x)
#B = np.array([np.cos(k*theta) for k in range(n+1)])

# direct evaluation of Legendre polynomials
B = legvander(x,n).T

# wall time elapsed
t2 = now()
print(f'Evaluation takes {t2-t1}s to complete')

#######################################################################

import matplotlib.pyplot as plt

plt.figure(figsize=(8,4))

plt.plot(x, P.T, '-', linewidth=3)
#plt.plot(x, P[n]-B[n], 'r-')

plt.axhline(0.0, color='black', linestyle=':', zorder=0)

plt.xlim([-1.0,1.0])

# show in interactive console
plt.show()
