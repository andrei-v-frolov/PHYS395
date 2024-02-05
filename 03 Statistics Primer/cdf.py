#!/usr/bin/env python
# rank IID random samples to estimate CDF and PDF

#######################################################################

# operating system functions (for standard stream access)
from sys import stdin

# numerical libraries
import numpy as np

#######################################################################

# read in data from stdin
data = np.loadtxt(stdin)

# sanity check on data format
assert len(data.shape) == 1, "Expecting single column data..."

# sorted data estimates the CDF
x = np.sort(data); rank = np.linspace(0.0, 1.0, len(x))

#######################################################################

# number of percentile bins
n = 32

# decimate the CDF using linear interpolation
theta = np.linspace(np.pi, 0.0, n)
cdf = (np.cos(theta) + 1.0)/2.0
x = np.interp(cdf, rank, x)

# compute the PDF from decimated CDF derivative
pdf = 1.0/np.gradient(x, cdf, edge_order=2)

#######################################################################

import matplotlib.pyplot as plt

#plt.plot(x, cdf, "r-")
plt.fill_between(x, pdf)

plt.ylim(bottom=0.0)

plt.show()
