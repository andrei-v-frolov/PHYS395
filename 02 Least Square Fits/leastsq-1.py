#!/usr/bin/env python
# generalized least square fit to supplied data
# run as: python leastsq.py < DATA
# expected data format: x f [sigma]

#######################################################################

# operating system functions (for standard stream access)
import sys

# numerical libraries
import numpy as np
from numpy.linalg import solve

#######################################################################

# read in data from stdin
data = np.loadtxt(sys.stdin)

# sanity check on supplied data format
n,columns = data.shape
assert columns >= 2, ("Expecting at least 2 columns, got %i" % columns)

# parsed data to be fitted
x = data[:,0]; f = data[:,1]

# optional weighting column
w = 1.0/data[:,2]**2 if columns > 2 else np.ones(n)

# basis functions used to fit the data
def basis(x):
	return np.array([1,x*x,x**4])

# initialize least square fit accumulators
A = np.zeros([3,3]); y = np.zeros(3)

# accumulate least square fit matrices
for i in range(0,n):
	b = basis(data[i,0])
	y += w[i]*f[i]*b
	A += w[i]*np.outer(b,b)

# solve for best fit coefficients
c = solve(A,y)

#######################################################################

import matplotlib.pyplot as plt

# compute best fit
y = np.array([np.dot(c,basis(x[i])) for i in range(0,n)])

# plot residual
plt.plot(x,f,'b-')
plt.plot(x,y,'r-')

plt.show()