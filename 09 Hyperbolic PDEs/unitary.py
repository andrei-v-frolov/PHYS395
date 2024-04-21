#!/usr/bin/env python
# unitary evolution of Schrodinger equation using finite differences

#######################################################################

import numpy as np
from scipy.signal import convolve
from scipy.linalg import solve_banded

#######################################################################

# grid resolution and box size
n = 1000; l = 5.0

# uniform spatial grid
dx = 2.0*l/(n-1); x = np.linspace(-l,l,n)

# initial wavefunction
psi = np.exp(-(x-1.5)**2/2.0).astype(complex)

#######################################################################

# time step
dt = dx

# oscillator potential
V = x*x/2.0

# kinetic term stencil
K = -np.array([1,-2,1])/(2.0*dx**2)

# band-diagonal Hamiltonian operator
H = np.repeat(K.reshape((3,1)), n, axis=1); H[1] += V

# implicit time step evolution matrix
Q = H*(0.5j*dt); Q[1] += 1.0

# advance solution to the next step (enforcing unitarity via Cayley transform)
def step(i):
	global psi
	# explicit half-step
	psi -= (convolve(psi, K, mode='same', method='direct') + V*psi)*(0.5j*dt)
	# implicit half-step
	psi = solve_banded((1,1), Q, psi)

#######################################################################

import matplotlib.pyplot as plt

fig = plt.figure(); ax = fig.gca()

p = plt.fill_between(x, np.abs(psi)**2, color='grey', alpha=0.5, linewidth=0.0)
re, = plt.plot(x, np.real(psi), 'r')
im, = plt.plot(x, np.imag(psi), 'b')

plt.xlim([-l,l])
plt.ylim([-1.1,1.1])

#######################################################################

import matplotlib.animation as animation

# called to advance animation to next frame
def animate(i):
	step(i+1)
	vertices = p.get_paths()[0].vertices
	vertices[1:n+1,1] = np.abs(psi)**2
	re.set_data(x, np.real(psi))
	im.set_data(x, np.imag(psi))

animation = animation.FuncAnimation(fig, animate, frames=1000, interval=1000.0/60)
#animation.save('unitary.mp4')

plt.show()
