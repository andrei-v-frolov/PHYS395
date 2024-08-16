#!/usr/bin/env python
# compute the covariance of a random process, and draw Gaussian sample with the same...
# run as: python mock.py < DATA

#######################################################################

# operating system functions (for standard stream access)
import sys

# numerical libraries
import numpy as np
from numpy.random import normal
from numpy.linalg import svd, eigh, cholesky

#######################################################################

# read in data from stdin
data = np.loadtxt(sys.stdin)

# sanity check on supplied data format
n,columns = data.shape
assert columns >= 2, ("Expecting at least 2 columns, got %i" % columns)

# normalize 1-point PDF (if required)
data = np.log(data)

# average value estimator
avg = np.sum(data,1)/columns

# initialize covariance accumulator
C = np.zeros([n,n])

# accumulate covariance matrix
for i in range(0,columns):
	x = data[:,i] - avg
	C += np.outer(x,x)

# normalize to # of samples
C /= columns-1

# compute correlation matrix
#D = np.diag(C); Q = C/np.sqrt(np.outer(D,D))

#######################################################################

# compute Cholesky decomposition C = L @ L.T
# (croaks for matrices of less than maximal rank)
#L = cholesky(C)

# compute eigenvalue decomposition C = U @ np.diag(S) @ U.T
#S,U = eigh(C); L = U * np.sqrt(np.maximum(0,S))

# compute SVD of the covariance C = U @ np.diag(S) @ V.T
U,S,VT = svd(C, hermitian=True); L = U * np.sqrt(S)

#######################################################################

# draw Gaussian samples with the same covariance
samples = 1500; X = np.dot(L, normal(size=[n,samples]))

#######################################################################

import matplotlib.pyplot as plt

#plt.imshow(C, origin='lower', cmap='inferno')
#plt.imshow(C, origin='lower', norm='log', cmap='inferno')
#plt.colorbar()

for i in range(0,samples):
	plt.plot(X[:,i]+avg, "r-", alpha=0.03)

for i in range(0,columns):
	plt.plot(data[:,i], "k-", alpha=0.50)

plt.xlim([0,n])
plt.show()
