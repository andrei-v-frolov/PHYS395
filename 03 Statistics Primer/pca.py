#!/usr/bin/env python
# compute SVD of supplied data, and keep only a few components
# run as: python pca.py < DATA

#######################################################################

# operating system functions (for standard stream access)
import sys

# numerical libraries
import numpy as np
from numpy.linalg import svd

#######################################################################

# read in data from stdin
data = np.loadtxt(sys.stdin)

# sanity check on supplied data format
n,columns = data.shape
assert columns >= 2, ("Expecting at least 2 columns, got %i" % columns)

#######################################################################

# compute SVD of the data
U,S,VT = svd(data, full_matrices=False)

# junk all but the most prominent components
S[5:] = 0.0

# reconstruct the data
X = U @ np.diag(S) @ VT

'''
# output it to stdout
for i in range(0,n):
	print(*X[i,:])
'''

#######################################################################

import matplotlib.pyplot as plt

plt.plot(X)
plt.legend(["AAPL", "AMD", "AMZN", "CSCO", "GOOGL", "INTC", "META", "MSFT", "NFLX", "NVDA", "QCOM", "SBUX", "TSLA"])
plt.xlim([0,2516])
plt.ylim([0,750])

#plt.scatter(X[:,0], X[:,11], c=np.linspace(0.0, 1.0, n), cmap='plasma', marker=".")
#plt.colorbar()

plt.show()
