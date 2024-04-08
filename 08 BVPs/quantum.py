#!/usr/bin/env python
# energy eigenstates of a quantum anharmonic oscillator

#######################################################################

from math import *
import numpy as np
from numpy.linalg import eigh
from numpy.polynomial.hermite import hermvander as Hn

#######################################################################

# basis size and highest eigenmode to plot
n = 100; m = n//8; xmax = 2.25*(m+1)**0.25

# Hamiltonian of a harmonic oscillator H = p^2/2 + x^2/2
H0 = np.diag(np.arange(n) + 0.5)

# position and momentum operator squared
X2 = H0.copy(); P2 = H0.copy()

for k in range(n-2):
	u = np.sqrt((k+1)*(k+2))/2.0
	X2[k,k+2] = u; P2[k,k+2] = -u
	X2[k+2,k] = u; P2[k+2,k] = -u

# spatial grid to evaluate wavefunctions on
x = np.linspace(-xmax, xmax, 1000)

# basis is eigenstates of a harmonic oscillator
B = Hn(x,n-1); w = np.exp(-x*x/2.0)/(np.pi**0.25)

for i in range(n):
	B[:,i] *= w/np.sqrt(2.0**i * factorial(i))

###############################################################################

# anharmonic oscillator potential
mu = 0.0; nu = 0.5; V = mu*x*x/2.0 + nu*x**4/4.0

# Hamiltonian H = p^2/2 + mu*x^2/2 + nu*x^4/4
H = P2/2.0 + mu*X2/2.0 + nu*np.matmul(X2,X2)/4.0

# compute the eigenvalues and eigenvectors of H
E,W = eigh(H)

# energy eigenstates of an anharmonic oscillator
psi = np.matmul(B,W)

###############################################################################

import matplotlib.pyplot as plt

# plot eigenstates and eigenvalues side-to-side
fig, ax = plt.subplots(1, 2, sharey=True); state = []

# plot PDF of energy eigenstates
for k in range(m):
	plot, = ax[0].plot(x, E[k]+psi[:,k]**2, color='black'); state.append(plot)

# plot oscillator potential
potential, = ax[0].plot(x, V, color='red', linewidth=3)
ax[0].set_xlim([-xmax,xmax])

# plot energy eigenvalues
energy, = ax[1].plot(range(m), E[:m], "o-")
ax[1].set_xlim([0,m-1])

plt.tight_layout()

#######################################################################

# interactive controls
from matplotlib.widgets import Slider

# make room for widgets
fig.subplots_adjust(bottom=0.22)

# user interface elements
mu_slider = Slider(
    ax=fig.add_axes([0.1, 0.07, 0.8, 0.03]),
    valmin=-2.0, valmax=2.0, valinit=mu,
    label='μ'
)

nu_slider = Slider(
    ax=fig.add_axes([0.1, 0.02, 0.8, 0.03]),
    valmin=0.0, valmax=1.0, valinit=nu,
    label='λ'
)

# update simulation parameters
def update(value):
	mu = mu_slider.val; nu = nu_slider.val
	H = P2/2.0 + mu*X2/2.0 + nu*np.matmul(X2,X2)/4.0
	E,W = eigh(H); psi = np.matmul(B,W)
	V = mu*x*x/2.0 + nu*x**4/4.0
	energy.set_data(range(m), E[:m])
	potential.set_data(x, V)
	for k in range(m):
		state[k].set_data(x, E[k]+psi[:,k]**2)
	ax[0].set_ylim([min(np.min(V), E[0])-0.5, E[m]])

# register update handler
mu_slider.on_changed(update)
nu_slider.on_changed(update)

update(None)
plt.show()
