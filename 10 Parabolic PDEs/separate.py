#!/usr/bin/env python
# phase separation modelled by Cahn-Hilliard equation

#######################################################################

import numpy as np
from scipy.ndimage import convolve
from numpy.random import choice

#######################################################################

# grid resolution and box size
n = 100; l = 1.0

# uniform spatial grid
dx = 2.0*l/n; x = np.linspace(-l,l,n)

# initial field profile
phi = choice([-1.0,1.0], (n,n))

#######################################################################

# second order Laplacian stencil (nearest neighbour)
#laplacian = np.array([[0,1,0],[1,-4,1],[0,1,0]])/dx**2

# second order Laplacian stencil (isotropic)
laplacian = np.array([[1,4,1],[4,-20,4],[1,4,1]])/(6.0*dx**2)

#######################################################################

# time step
dt = dx**4

# advance solution to the next step (enforcing BCs)
def step(i):
	global phi
	mu = phi**3 - phi - convolve(phi, laplacian, mode='wrap') * dx**2
	phi += convolve(mu, laplacian, mode='wrap') * dt
	return phi

#######################################################################

import matplotlib.pyplot as plt
from matplotlib import colors

fig = plt.figure(); ax = fig.gca()
crt = colors.LinearSegmentedColormap.from_list("CRT", ['black', 'greenyellow'])
wave = plt.imshow(phi, extent=[-l,l,-l,l], vmin=-1.0, vmax=1.0, cmap=crt, interpolation='none')

#######################################################################
import matplotlib.animation as animation

# called to advance animation to next frame
def animate(i):
	for k in range(0,n): step(i)
	wave.set_data(phi)
	#z = np.max(phi); wave.set_clim([-z,z])

animation = animation.FuncAnimation(fig, animate, frames=1000, interval=1000.0/60)
#animation.save('separate.mp4')

plt.show()
