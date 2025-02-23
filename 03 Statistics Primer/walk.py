#!/usr/bin/env python
# stochastic process (aka random walk) demo

#######################################################################

import numpy as np
from numpy.random import normal

# number of steps and realizations
n = 1024; chains = 500

# draw all the random numbers we need at once
step = normal(size=[n,chains])

# fast way to sum random walk steps
#x = np.cumsum(step, axis=0)

# otherwise, we can do it ourselves...
x = np.zeros_like(step); x[0] = 1.0

for i in range(0,n-1):
	x[i+1] = (1.002 + 0.015*step[i])*x[i]

#######################################################################

import matplotlib.pyplot as plt

plt.plot(x, 'r-', alpha=0.03)
plt.plot(x[:,0], 'b-')
plt.xlim([0,n])

plt.show()
