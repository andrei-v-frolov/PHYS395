#!/usr/bin/env python
# heat diffusion equation (in two dimensions)

#######################################################################

import numpy as np
from scipy.ndimage import convolve

#######################################################################

# grid resolution and box size
n = 100; l = 1.0

# uniform spatial grid
dx = 2.0*l/(n-1); x = np.linspace(-l,l,n)

# initial field profile
X,Y = np.meshgrid(x,x)
phi = 3.0*np.exp(-(X*X+Y*Y)*256.0)

#######################################################################

# second order Laplacian stencil (nearest neighbour)
#laplacian = np.array([[0,1,0],[1,-4,1],[0,1,0]])/dx**2

# second order Laplacian stencil (isotropic)
laplacian = np.array([[1,4,1],[4,-20,4],[1,4,1]])/(6.0*dx**2)

#######################################################################

# time step
dt = dx**2/4.0

# boundary conditions
def bc():
	phi[0,:] = 1.0; phi[-1,:] = 0.0
	#phi[:,0] = 1.0; phi[:,-1] = 0.0

# pre-computed wave equation stencil
stencil = np.array([[0,0,0],[0,1,0],[0,0,0]]) + laplacian*dt

# advance solution to the next step (enforcing BCs)
def step(i):
	global phi
	phi = convolve(phi, stencil, mode='constant'); bc()
	return phi

#######################################################################

import matplotlib.pyplot as plt

fig = plt.figure(); ax = fig.gca()
wave = plt.imshow(phi, extent=[-l,l,-l,l], vmin=-1.0, vmax=1.0, cmap='seismic')

#######################################################################
import matplotlib.animation as animation

# called to advance animation to next frame
def animate(i):
	wave.set_data(step(i))
	#z = np.max(phi); wave.set_clim([-z,z])

animation = animation.FuncAnimation(fig, animate, frames=1000, interval=1000.0/60)
#animation.save('diffuse.mp4')

plt.show()