#!/usr/bin/env python
# compute the covariance of supplied data, and draw Gaussian sample with the same...
# run as: python covariance.py < DATA > GAUSSIAN

#######################################################################

# operating system functions (for standard stream access)
import sys

# numerical libraries
import numpy as np
from numpy.random import normal
from numpy.linalg import svd, cholesky

#######################################################################

# read in data from stdin
data = np.loadtxt(sys.stdin)

# sanity check on supplied data format
n,columns = data.shape
assert columns >= 2, ("Expecting at least 2 columns, got %i" % columns)

# subtract average from the data (aka de-mean)
avg = np.sum(data,0)/n; data -= avg

# initialize covariance accumulator
C = np.zeros([columns,columns])

# accumulate covariance matrix
for i in range(0,n):
	x = data[i,:]
	C += np.outer(x,x)

# normalize to # of samples
C /= n

# compute correlation matrix
D = np.diag(C); Q = C/np.sqrt(np.outer(D,D))

#######################################################################

# draw Gaussian samples with the same covariance
L = cholesky(C); X = normal(size=[n,columns]) @ L.T

# output them to stdout
for i in range(0,n):
	print(*(X[i,:]+avg))

#######################################################################

'''
# compute SVD of the covariance matrix
U,S,VT = svd(C, hermitian=True)

# junk all but the most prominent component
S[1:] = 0.0

# reconstruct the covariance matrix
C = (U*S) @ VT
'''

#######################################################################

import matplotlib.pyplot as plt

plt.imshow(C, origin='lower', cmap='inferno')
#plt.imshow(C, origin='lower', norm='log', cmap='inferno')

plt.colorbar()
plt.show()
