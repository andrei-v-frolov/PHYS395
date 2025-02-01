#!/usr/bin/env python
# generate some test data for least square fit
# run as: python model.py > DATA

#######################################################################

import numpy as np
from numpy.random import normal, uniform

#######################################################################

# number of points to generate
n = 371

# uniform or random sampling
#x = np.linspace(-1.0,1.0,n)
x = np.sort(uniform(-1.0,1.0,n))

# test function
f = np.exp(-x*x*4.5)

# add noise to make things more interesting
f += normal(scale=0.1, size=n)

# output whitespace-separated formatted data
for i in range(0,n):
	print("%24.16f %24.16f" % (x[i], f[i]))
