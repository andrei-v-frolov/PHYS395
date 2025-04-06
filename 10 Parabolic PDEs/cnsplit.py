#!/usr/bin/env python
# heat diffusion equation using Crank-Nicolson method

#######################################################################

import numpy as np
from scipy.signal import convolve
from scipy.linalg import solve_banded

#######################################################################

# grid resolution and box size
n = 256; l = 1.0

# uniform spatial grid
dx = 2.0*l/(n-1); x = np.linspace(-l,l,n)

# initial field profile
X,Y = np.meshgrid(x,x)
phi = 3.0*np.exp(-(X*X+Y*Y)*256.0)

#######################################################################

# second order Laplacian stencil in one dimension
laplacian = np.array([1,-2,1])/dx**2

#######################################################################

# time step
dt = dx/4.0

# pre-computed forward and backward half-step stencils
fwd = np.array([0,1,0]) + laplacian*(dt/2.0)
bwd = np.array([0,1,0]) - laplacian*(dt/2.0)

# implicit half-step solver matrix
Q = np.repeat(bwd.reshape((3,1)), n, axis=1)

# advance solution to the next step (enforcing BCs)
def step(i):
	for i in range(n): phi[i,:] = convolve(phi[i,:], fwd, mode='same', method='direct')
	for i in range(n): phi[:,i] = convolve(solve_banded((1,1), Q, phi[:,i]), fwd, mode='same', method='direct')
	for i in range(n): phi[i,:] = solve_banded((1,1), Q, phi[i,:])
	return phi

#######################################################################

import matplotlib.pyplot as plt

fig = plt.figure(); ax = fig.gca()
wave = plt.imshow(phi, cmap='seismic', extent=[-l,l,-l,l], origin='lower', interpolation='none')
plt.colorbar()

#######################################################################

import matplotlib.animation as animation

# called to advance animation to next frame
def animate(i):
	wave.set_data(step(i))
	z = np.max(phi); wave.set_clim([-z,z])

animation = animation.FuncAnimation(fig, animate, frames=1000, interval=1000.0/60)
#animation.save('diffuse.mp4')

plt.show()
