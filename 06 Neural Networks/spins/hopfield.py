#!/usr/bin/env python
# Hopfield (aka associative memory) spin network

#######################################################################

import numpy as np
from numpy.random import rand, choice

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
img = Image.open('digits.png')
alpha = np.array(img)[...,-1]/255.0

# slice and threshold artboard into individual digits
digit = [np.where(alpha[:,28*i:28*(i+1)] > 0.75, 1.0, -1.0) for i in range(10)]

#######################################################################

# spin coupling matrix
W = np.zeros([n*n,n*n])

# Hebbian learning rule
for i in range(10):
	W += np.outer(digit[i],digit[i])

# remove self-coupling
np.fill_diagonal(W, 0.0)

# synchronous state update
def update(sigma):
	return np.where(np.dot(W,sigma.flat).reshape(n,n) > theta, 1.0, -1.0)

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
