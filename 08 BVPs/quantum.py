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

# physical wavefunction normalization
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

# canonical sign choice for eigenstates
for i in range(n):
	psi[:,i] *= np.sign(np.dot(1.0+x,psi[:,i]))

###############################################################################

import matplotlib.pyplot as plt

# plot eigenstates and eigenvalues side-by-side
fig, ax = plt.subplots(1, 2, sharey=True)

# what to plot for energy level
plot = 'probability'; phase = np.zeros(m); state = []

def level(k):
	return (psi[:,k]**2 if (plot == 'probability') else psi[:,k]*cos(phase[k])/2.0) + E[k]

# plot energy eigenstates
for k in range(m):
	lvl, = ax[0].plot(x, level(k), color='black'); state.append(lvl)

# plot oscillator potential
potential, = ax[0].plot(x, V, color='red', linewidth=3)
ax[0].set_xlim([-xmax,xmax])

# plot energy eigenvalues
energy, = ax[1].plot(range(m), E[:m], "o-")
ax[1].set_xlim([0,m-1])

plt.tight_layout()

#######################################################################

# interactive controls
from matplotlib.widgets import Slider, RadioButtons

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

choice = RadioButtons(
	fig.add_axes([0.65, 0.83, 0.22, 0.12]),
	('probability','wavefunction','animation')
)

# change display mode
def change(value):
	global plot; plot = choice.value_selected
	if (plot != 'animation'): phase.fill(0.0)
	if (plot == 'wavefunction'):
		# reset to canonical sign
		psi = np.expand_dims(1.0+x,1)
	update(value)

# update simulation parameters
def update(value):
	global mu, nu, V, H, E, W, psi
	# update model parameters
	mu = mu_slider.val; nu = nu_slider.val
	H = P2/2.0 + mu*X2/2.0 + nu*np.matmul(X2,X2)/4.0
	V = mu*x*x/2.0 + nu*x**4/4.0
	# recompute eigenfunctions inheriting sign
	E,W = eigh(H); phi = np.matmul(B,W)
	for i in range(n):
		psi[:,i] = np.sign(np.dot(psi[:,i],phi[:,i]))*phi[:,i]
	# update plot data
	energy.set_data(range(m), E[:m])
	potential.set_data(x, V)
	for k in range(m): state[k].set_data(x, level(k))
	ax[0].set_ylim([min(np.min(V), E[0])-0.5, E[m]])

# register update handler
mu_slider.on_changed(update)
nu_slider.on_changed(update)
choice.on_clicked(change)

update(None)

#######################################################################

import matplotlib.animation as animation

# evolution time step
dt = 4.0*np.pi/1000

# called to advance animation to next frame
def animate(i):
	if (plot != 'animation'): return
	global phase; phase += E[:m]*dt
	for k in range(m):
		state[k].set_data(x, level(k))

animation = animation.FuncAnimation(fig, animate, frames=1000, interval=1000.0/60)
#animation.save('quantum.mp4')
plt.show()
