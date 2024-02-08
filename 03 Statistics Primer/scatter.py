#!/usr/bin/env python
# scatter plot of bivariate PDF, with a few variations

#######################################################################

# operating system functions (for standard stream access)
from sys import stdin

# numerical libraries
import numpy as np
from erfinv import erfinv
from scipy.stats import rankdata
from numpy.linalg import svd

#######################################################################

# read in data from stdin
data = np.loadtxt(stdin)
n, columns = data.shape

# sanity check on data format
assert columns > 1, "Expecting at least two columns of data..."

#######################################################################

import matplotlib.pyplot as plt

x = data[:,0]	# Apple
y = data[:,11]	# Starbucks
z = data[:,8]	# Netflix

t = np.linspace(0.0, 1.0, n)

#######################################################################

'''
# equalize (or normalize) the univariate distributions
x = (rankdata(x)-0.5)/n; x = np.sqrt(2.0) * np.vectorize(erfinv)(2.0*x-1.0)
y = (rankdata(y)-0.5)/n; y = np.sqrt(2.0) * np.vectorize(erfinv)(2.0*y-1.0)
'''

#######################################################################

'''
# de-mean the data
x -= np.sum(x)/n
y -= np.sum(y)/n

# covariance components
xx = np.sum(x*x)/n
xy = np.sum(x*y)/n
yy = np.sum(y*y)/n

# correlation coefficient
c = xy/np.sqrt(xx*yy)

# SVD decomposition of covariance (symmetric, hence VT = U)
U,S,VT = svd(np.array([[xx,xy],[xy,yy]]), hermitian=True)

# whiten the data using covariance decomposition
X = np.row_stack((x,y))
X = np.matmul(U.T,X)
X = X/np.sqrt(S).reshape(2,1)
X = np.matmul(U,X)

x = X[0,:]
y = X[1,:]
'''

#######################################################################

# scatter plot, with (optional) color coding
#plt.scatter(x, y, c=t, cmap='plasma', marker=".")
plt.scatter(x, y, c=z, cmap='plasma', marker=".")

# binned histograms, with log(pdf) used for color
#plt.hist2d(x, y, bins=30, norm='log', cmap='Blues')
#plt.hexbin(x, y, gridsize=30, bins='log', cmap='Blues')

# set square aspect if desired
#plt.gca().set_aspect('equal')

plt.colorbar()
plt.show()
