#!/usr/bin/env python
# Ising model in two dimensions

#######################################################################

import numpy as np
from scipy.ndimage import convolve
from numpy.random import rand, choice

#######################################################################

# grid resolution
n = 101

# site parity mask
X,Y = np.meshgrid(range(n),range(n)); parity = X+Y

# magnetic field template
field = np.ones([n,n])

# initial spin configuration
sigma = choice([-1.0,1.0], (n,n))

#######################################################################

# coupling constants
J = 1.0; mu = 0.0; beta = 1.0

# nearest neighbour coupling stencil
stencil = np.array([[0,1,0],[1,0,1],[0,1,0]])

#######################################################################

# split Metropolis-Hastings MCMC update
def step(i):
	global sigma
	force = J*convolve(sigma, stencil, mode='wrap') + mu*field
	trial = choice([-1.0,1.0], (n,n)); alpha = np.exp(beta*force*(trial-sigma))
	sigma = np.where((rand(n,n) > alpha) | ((parity+i)%2 == 0), sigma, trial)
	return sigma

#######################################################################

import matplotlib.pyplot as plt
from matplotlib import colors

crt = colors.LinearSegmentedColormap.from_list("CRT", ["black", "greenyellow"])

fig = plt.figure(); ax = fig.gca()
spins = plt.imshow(sigma, origin='lower', vmin=-1.0, vmax=1.0, cmap=crt)

#######################################################################
import matplotlib.animation as animation

# called to advance animation to next frame
def animate(i):
	spins.set_data(step(i))

animation = animation.FuncAnimation(fig, animate, frames=1000, interval=1000.0/60)
#animation.save('ising.mp4')

plt.show()