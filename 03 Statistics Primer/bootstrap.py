#!/usr/bin/env python
# re-draw samples (drawn from IID) randomly using bootstrap

#######################################################################

# operating system functions (for standard stream access)
from sys import stdin

# numerical libraries
import numpy as np
from numpy.random import rand

#######################################################################

# read in data from stdin
data = np.loadtxt(stdin)

# sanity check on data format
assert len(data.shape) == 1, "Expecting single column data..."

# random indices selecting sub-sample of data
n = len(data); idx = np.floor(n*rand(n))

# output samples, one per line, in scientific notation
for i in range(0,n):
	print("%.16e" % data[int(idx[i])])
