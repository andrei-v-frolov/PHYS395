#!/usr/bin/env python
# Hopfield (aka associative memory) spin network

#######################################################################

import numpy as np
from numpy.random import rand, choice

# default random number generator
rng = np.random.default_rng()

#######################################################################

# grid resolution
n = 28

# external field template
theta = np.zeros([n,n])

# initial spin configuration
sigma = choice([-1.0,1.0], (n,n))

#######################################################################

from PIL import Image

# import 280x28 artboard with target digit images
img = np.array(Image.open('digits.png'))

# slice and threshold artboard into individual digits
digit = [np.where(img[:,i*n:(i+1)*n,-1] > 190, 1.0, -1.0) for i in range(10)]

#######################################################################

# flatten 2D images into state vectors
V = np.array(digit).reshape(10,n*n)

# Hebbian learning rule (does not work that well...)
#W = V.T @ V

# use de-correlated vectors instead (kind of works)
C = V @ V.T; W = V.T @ np.linalg.inv(C) @ V

# remove self-coupling (not that it matters for unit spins)
np.fill_diagonal(W, 0.0)

#######################################################################

# synchronous downhill state update
def update(sigma):
	return np.where(np.dot(W,sigma.flat).reshape(n,n) > theta, 1.0, -1.0)

# synchronous MCMC state update
def mcmc(sigma, beta=2.0):
	s = sigma.flat
	t = choice([-1.0,1.0], n*n)
	force = np.dot(W,s) - theta.flat
	alpha = np.exp(beta*force*(t-s))
	return np.where(rand(n*n) > alpha, s, t).reshape(n,n)

#######################################################################

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap as colormap

fig = plt.figure(); ax = fig.gca()

# spin orientation map
crt = colormap.from_list("CRT", ['black', 'greenyellow'])
spins = plt.imshow(sigma, origin='upper', cmap=crt, vmin=-1.0, vmax=1.0, interpolation='none')

# no ticks and tight framing
plt.tick_params(left=False, right=False, labelleft=False, labelbottom = False, bottom=False)
plt.tight_layout()

#######################################################################

import matplotlib.animation as animation

# called to advance animation to next frame
def animate(i):
	global sigma
	sigma = update(sigma)
	spins.set_data(sigma)

animation = animation.FuncAnimation(fig, animate, frames=1000, interval=1000.0/60)
#animation.save('hopfield.mp4')

plt.show()
