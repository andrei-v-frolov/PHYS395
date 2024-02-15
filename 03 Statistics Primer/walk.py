#!/usr/bin/env python
# stochastic processes (aka random walk) demo

#######################################################################

# number of steps and realizations
n = 1024; chains = 500

import numpy as np
from numpy.random import normal, uniform

# draw all the random numbers we need at once
step = normal(size=[n,chains])

# fast way to sum random walk steps
#x = np.cumsum(step, axis=0)

# otherwise, we can do it ourselves...
x = np.zeros([n,chains]); x[:,:] = 1.0

for i in range(0,n-1):
	x[i+1] = x[i] + 0.002*x[i] + 0.015*x[i]*step[i]

#######################################################################

import matplotlib.pyplot as plt

plt.plot(x, "r-", alpha=0.03)
plt.plot(x[:,0], "b-")
plt.show()
