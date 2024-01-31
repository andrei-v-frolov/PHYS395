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
from numpy.polynomial.legendre import legvander, legval

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
B = legvander(x,4)
A = np.matmul(B.transpose(),B)
y = np.matmul(B.transpose(),f)

# solve bor best fit coefficients
c = solve(A,y)

#######################################################################

import matplotlib.pyplot as plt

# compute best fit
y = legval(x,c)

# plot residual
plt.plot(x,f,'b-')
plt.plot(x,y,'r-')

plt.show()