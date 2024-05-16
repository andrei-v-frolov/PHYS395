#!/usr/bin/env python
# MCMC sampler using Metropolis-Hastings algorithm

#######################################################################

import numpy as np
from numpy.random import normal, uniform

#######################################################################

# number of steps and realizations
n = 1024; chains = 1000; l = 1.5

# initial chain positions
x = uniform(-l, l, size=[2,chains])

# slightly tilted mexican hat likelihood
def likelihood(x,y):
	U = 20.0*(x*x + y*y - 1.0)**2 + y
	return np.exp(-U)

'''
# random Hamiltonian kick, with dt ~ T/32
def kick(state, pi=5.0, dt=0.015, steps=16):
	x,y = state.copy()
	vx,vy = normal(scale=pi, size=[2,chains])
	for k in range(steps):
		x += vx * (dt/2.0); y += vy * (dt/2.0)
		K = 80.0*(x*x + y*y - 1.0)
		vx -= (K*x) * dt; vy -= (K*y+1.0) * dt
		x += vx * (dt/2.0); y += vy * (dt/2.0)
	return np.array([x,y])
'''

# random step size (adjust for acceptance rate ~ 0.5)
sigma = 0.3

# MCMC step using Metropolis-Hastings algorithm
def step(x):
	#y = kick(x)
	y = x + normal(scale=sigma, size=[2,chains])
	alpha = likelihood(y[0], y[1])/likelihood(x[0], x[1])
	
	u = uniform(size=chains) > alpha
	#print('acceptance rate = %.3f' % (1.0 - np.sum(u)/chains))
	return np.where(np.reshape(u, [1,chains]), x, y)

#for i in range(0,n): x = step(x)

#######################################################################

import matplotlib.pyplot as plt

# sample likelihood on a grid
grid = np.linspace(-l, l, 1024)
X,Y = np.meshgrid(grid, grid)
L = likelihood(X,Y)

fig = plt.figure(); ax = fig.gca()

# density plot of likelihood
plt.imshow(L, cmap='Blues', origin='lower', extent=[-l,l,-l,l], interpolation='none')
#plt.contourf(X,Y,L, levels=30, cmap='Blues')
plt.colorbar()

# current state of the chains
walkers = plt.scatter(x[0],x[1], marker=".", color="tab:red", s=5.0)
#walkers = plt.scatter(x[0],x[1], c=np.arctan2(x[0],x[1]), cmap='plasma', marker=".", s=5.0)

# plot extent
ax.set_aspect('equal')
plt.xlim([-l,l])
plt.ylim([-l,l])

#######################################################################

import matplotlib.animation as animation

# called to advance animation to next frame
def animate(i):
	global x; x = step(x)
	walkers.set_offsets(x.T)

animation = animation.FuncAnimation(fig, animate, frames=1000, interval=1000/60)
#animation.save('mcmc.mp4')
plt.show()
