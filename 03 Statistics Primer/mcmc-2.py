#!/usr/bin/env python
# MCMC sampler using Metropolis-Hastings algorithm

#######################################################################

# number of steps and realizations
n = 1024; chains = 500; l = 1.5

import numpy as np
from numpy.random import normal, uniform

# initial chain positions
x = uniform(-l, l, size=[2,chains])

# slightly tilted mexican hat likelihood
def likelihood(x,y):
	v = (x*x + y*y - 1.0)**2 + 0.05*y
	return np.exp(-20.0*v)

# random step size (adjust for acceptance rate ~ 0.5)
sigma = 0.3

# store chain history (you might want to downsample)
history = np.zeros([n,2,chains])

# MCMC using Metropolis-Hastings algorithm
for i in range(0,n):
	y = x + normal(scale=sigma, size=[2,chains])
	alpha = likelihood(y[0], y[1])/likelihood(x[0], x[1])
	
	u = uniform(size=chains) > alpha
	x = np.where(np.reshape(u, [1,chains]), x, y)
	history[i] = x

#######################################################################

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation as animate

# sample likelihood on a grid
grid = np.linspace(-l, l, 1024)
X,Y = np.meshgrid(grid, grid)
L = likelihood(X,Y)

fig = plt.figure()

# density plot of likelihood
plt.contourf(X,Y,L, levels=30, cmap='Reds')
plt.colorbar()

# end state of the chains
walkers = plt.scatter(x[0],x[1], marker=".")

#animation = animate(fig, lambda i: walkers.set_offsets(history[i].T), n, interval=1000/60)

plt.xlim([-l,l]); plt.ylim([-l,l])
plt.gca().set_aspect('equal')
#animation.save('mcmc.mp4')
plt.show()
