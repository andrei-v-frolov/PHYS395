#!/usr/bin/env python
# sort IID random samples to estimate CDF

#######################################################################

# standard input access
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

# decimate CDF using linear interpolation
cdf = np.sin(np.linspace(0.0, np.pi/2, n))**2
x = np.interp(cdf, rank, x)

# compute PDF from decimated CDF derivative
pdf = np.gradient(cdf, x)

#######################################################################

import matplotlib.pyplot as plt

#plt.plot(x, cdf)
plt.fill_between(x, pdf)

# restrict axis range
plt.xlim([x[0],x[-1]])
plt.ylim(bottom=0.0)

# show in interactive console
plt.show()
