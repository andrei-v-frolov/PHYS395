#/usr/bin/env python
# compute principal components of supplied data
# run as: python pca.py < DATA

#######################################################################

# standard stream access
from sys import stdin

# numerical libraries
import numpy as np
from numpy.linalg import svd, norm

#######################################################################

# read in data from stdin
data = np.loadtxt(stdin)
n,columns = data.shape

# sanity check on supplied data format
assert columns > 1, f'Expecting at least 2 columns, got {columns}'

#######################################################################

# compute singular value decomposition
U,S,V = svd(data, full_matrices=False)
print(f'SVD accuracy {norm(U*S@V - data)}')

# junk all but the most prominent components
S[7:] = 0.0

# reconstruct the data
model = (U*S)@V

'''
# output it to stdout
for i in range(0,n):
	print(*model[i,:])
'''

#######################################################################

import matplotlib.pyplot as plt

plt.plot(model)
plt.legend(["AAPL", "AMD", "AMZN", "CSCO", "GOOGL", "INTC", "META", "MSFT", "NFLX", "NVDA", "QCOM", "SBUX", "TSLA"], frameon=False)
plt.xlim([0,n])

#plt.scatter(model[:,0],model[:,11], c=np.arange(n), cmap='plasma', marker=".", s=3)
#plt.colorbar()

plt.show()
