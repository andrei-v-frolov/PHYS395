#!/usr/bin/env python
# visualize streamlines with line integral convolution

#######################################################################

import numpy as np
from numpy.random import normal
from scipy.ndimage import convolve

#######################################################################

# grid resolution and box size
n = 1024; l = 2.0

# uniform spatial grid
dx = 2.0*l/(n-1); x = np.linspace(-l,l,n)

# electric field of a dipole
X,Y = np.meshgrid(x,x)
X1 = X-1.0; X2 = X+1.0
R1 = (X1*X1+Y*Y)**1.5
R2 = (X2*X2+Y*Y)**1.5
Ex = X1/R1 - X2/R2
Ey = Y/R1 - Y/R2

#######################################################################

# anisotropic diffusion metric
E2 = Ex*Ex+Ey*Ey
uxx = Ex*Ex/E2
uxy = Ex*Ey/E2
uyy = Ey*Ey/E2

# time step
dt = dx**2/4.0

# pre-computed time step operators
Dxx = np.array([[1,-2,1]])*(dt/dx**2)
Dyy = np.array([[1],[-2],[1]])*(dt/dx**2)
Dxy = np.array([[1,0,-1],[0,0,0],[-1,0,1]])*(2.0*dt)/(4.0*dx**2)
L = np.array([[1,4,1],[4,-20,4],[1,4,1]])*(dt/4.0)/(6.0*dx**2)

#######################################################################

# initial seed and boundary conditions
phi = normal(size=(n,n)); bc = 'reflect'

# limit random field bandwidth to avoid artefacts
for i in range(16): phi += convolve(phi, L, mode=bc)

# advance solution to the next step (enforcing BCs)
def step(i):
	global phi
	phi += uxx*convolve(phi, Dxx, mode=bc) + uxy*convolve(phi, Dxy, mode=bc) + uyy*convolve(phi, Dyy, mode=bc)
	return phi

#######################################################################

import matplotlib.pyplot as plt

fig = plt.figure(); fig.gca().set_aspect('equal')
line = plt.imshow(phi, extent=[-l,l,-l,l], vmin=-1.0, vmax=1.0, cmap='gray', interpolation='none')
#plt.imshow(np.arcsinh(R1**(-1.0/3.0)-R2**(-1.0/3.0)), extent=[-l,l,-l,l], vmin=-3.0, vmax=3.0, cmap='bwr', interpolation='none')
#plt.streamplot(X,Y,Ex,Ey)

#######################################################################
import matplotlib.animation as animation

# called to advance animation to next frame
def animate(i):
	line.set_data(step(i))
	z = np.max(phi)/2.0; line.set_clim([-z,z])

animation = animation.FuncAnimation(fig, animate, frames=1000, interval=1000.0/60)
#animation.save('stream.mp4')

plt.show()
