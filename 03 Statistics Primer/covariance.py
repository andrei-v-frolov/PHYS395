#!/usr/bin/env python
# compute covariance of the supplied data, and (optionally)
#   - whiten the data (in several ways)
#   - draw Gaussian samples with the same...
# run as: python covariance.py < DATA [> OUTPUT]

#######################################################################

# standard stream access
from sys import stdin

# numerical libraries
import numpy as np
from numpy.random import normal
from numpy.linalg import svd, cholesky

#######################################################################

# read in data from stdin
data = np.loadtxt(stdin)
n,columns = data.shape

# sanity check on supplied data format
assert columns > 1, f'Expecting at least 2 columns, got {columns}'

#######################################################################

# subtract out the average (aka de-mean)
avg = np.sum(data,axis=0)/n; data -= avg

# covariance matrix
C = (data.T @ data)/(n-1)

# correlation matrix
D = np.sqrt(np.diag(C)); Q = C/np.outer(D,D)

#######################################################################

'''
# whiten data using covariance
U,S,V = svd(C, hermitian=True)
data = data @ (U/np.sqrt(S)@V)

# whiten data using correlation
#U,S,V = svd(Q, hermitian=True)
#data = data/D @ (U/np.sqrt(S)@V)

# whiten data using Cholesky decomposition
#U,S,V = svd(C, hermitian=True)
#data = data @ cholesky(U/S@V)

# output it to stdout
for i in range(n):
	print(*data[i])
'''

#######################################################################

'''
# draw Gaussian samples with the same covariance
L = cholesky(C); draw = normal(size=[n,columns]) @ L.T

# output them to stdout
for i in range(n):
	print(*(draw[i]+avg))
'''

#######################################################################

import matplotlib.pyplot as plt

plt.imshow(C, origin='lower', norm='log', cmap='inferno')
#plt.imshow(Q, origin='lower', cmap='inferno')

# scatter plot, with (optional) color coding
#plt.scatter(draw[:,0]+avg[0], draw[:,11]+avg[11], color='tab:blue', marker=".", s=3)
#plt.scatter(data[:,0]+avg[0], data[:,11]+avg[11], c=np.arange(n), cmap='plasma', marker=".", s=3)
#plt.scatter(data[:,0], data[:,11], c=np.arange(n), cmap='plasma', marker=".", s=3)

# set square aspect if desired
#plt.gca().set_aspect('equal')
#plt.xlim([-5,5]); plt.ylim([-5,5])

plt.colorbar()
plt.show()
