#!/usr/bin/env python
# wave equation using leapfrog method (in two dimensions)

#######################################################################

import numpy as np
from scipy.ndimage import convolve

#######################################################################

# grid resolution and box size
n = 1000; l = 1.0

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

# fancy boundary conditions (if desired)
#mask = (X*X + (X*X+2*Y)**4/4) > 0.9

#######################################################################

# time step
dt = dx/2.0

# ring buffer sampling field values
smp = np.zeros([3,n,n])

# momentarily stationary initial conditions
smp[0] = phi
smp[1] = phi

# pre-computed wave equation stencil
stencil = np.array([[0,0,0],[0,2,0],[0,0,0]]) + laplacian*dt**2

# advance solution to the next step (BCs encoded in mode)
def step(i):
	dn = i%3; hr = (i+1)%3; up = (i+2)%3
	smp[up] = convolve(smp[hr], stencil, mode='constant') - smp[dn]
	#smp[up] = np.where(mask, smp[hr], convolve(smp[hr], stencil, mode='constant') - smp[dn])
	return smp[up]

#######################################################################

import matplotlib.pyplot as plt

fig = plt.figure(); ax = fig.gca()
wave = plt.imshow(smp[1], cmap='seismic', vmin=-1.0, vmax=1.0, extent=[-l,l,-l,l], origin='lower', interpolation='none')
#plt.contourf(mask, extent=[-l,l,-l,l], colors=['#00000000', 'skyblue'])

#######################################################################

import matplotlib.animation as animation

# called to advance animation to next frame
def animate(i):
	wave.set_data(step(i))

animation = animation.FuncAnimation(fig, animate, frames=1000, interval=1000.0/60)
#animation.save('wave-2.mp4')

plt.show()