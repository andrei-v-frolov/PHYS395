#!/usr/bin/env python
# wave equation using leapfrog method (in one dimension)

#######################################################################

import numpy as np
from scipy.signal import convolve

#######################################################################

# grid resolution and box size
n = 1000; l = 1.0

# uniform spatial grid
dx = 2.0*l/(n-1); x = np.linspace(-l,l,n)

# initial field profile
phi = np.exp(-x*x*256.0)

#######################################################################

# time step
dt = dx

# ring buffer sampling field values
smp = np.zeros([3,n])

# momentarily stationary initial conditions
smp[0] = phi
smp[1] = phi

# pre-computed Laplacian evolution stencil
stencil = np.array([0,2,0]) + np.array([1,-2,1]) * (dt/dx)**2

# advance solution to the next step (enforcing specified BCs)
def step(i):
	dn = i%3; hr = (i+1)%3; up = (i+2)%3
	smp[up] = convolve(smp[hr], stencil, mode='same', method='direct') - smp[dn]
	smp[up,0] = 0.0; smp[up,-1] = smp[up,-2]
	return smp[up]

#######################################################################

import matplotlib.pyplot as plt

fig = plt.figure(); ax = fig.gca()
wave, = plt.plot(x, smp[1], 'r')

plt.xlim([-l,l])
plt.ylim([-1.1,1.1])

#######################################################################

import matplotlib.animation as animation

# called to advance animation to next frame
def animate(i):
	wave.set_data(x, step(i))

animation = animation.FuncAnimation(fig, animate, frames=1000, interval=1000.0/60)
#animation.save('wave-1.mp4')

plt.show()