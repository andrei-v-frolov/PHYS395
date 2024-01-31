#!/usr/bin/env python
# generate some test data for least square fit
# run as: python model.py > DATA

#######################################################################

import numpy as np

#######################################################################

# number of points to generate
n = 371

# test case - cosine on uniform grid
x = np.linspace(-1.0,1.0,n)
f = np.cos(np.pi*x)

# output whitespace-separated formatted data
for i in range(0,n):
	print("%24.16f %24.16f" % (x[i], f[i]))
