#!/usr/bin/env python
# MCMC sampler using Metropolis-Hastings algorithm

#######################################################################

# number of steps and realizations
n = 1024; chains = 500

import numpy as np
from numpy.random import normal, uniform

# initial chain positions
x = uniform(-15.0, 15.0, size=chains)

# sample Gaussian PDF (for starters)
def likelihood(x):
	return np.exp(-x*x/2.0)

# random step size (adjust for acceptance rate ~ 0.5)
sigma = 0.3

# store chain history (you might want to downsample)
history = np.zeros([n,chains])

for i in range(0,n):
	y = x + normal(scale=sigma, size=chains)
	alpha = likelihood(y)/likelihood(x)
	
	u = uniform(size=chains)
	x = np.where(u > alpha, x, y)
	history[i] = x

#######################################################################

import matplotlib.pyplot as plt

plt.plot(history, "r-", alpha=0.03)
plt.plot(history[:,0], "b-")
plt.xlim([0,n])
plt.show()
