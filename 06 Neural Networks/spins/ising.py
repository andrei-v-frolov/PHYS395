#!/usr/bin/env python
# Ising model in two dimensions

#######################################################################

import numpy as np
from scipy.ndimage import convolve
from numpy.random import rand, choice

#######################################################################

# grid resolution
n = 512

# site parity mask
X,Y = np.meshgrid(range(n),range(n)); parity = X+Y

# magnetic field template
field = np.ones([n,n])

# initial spin configuration
sigma = choice([-1.0,1.0], (n,n))

#######################################################################

# coupling constants
J = 1.0; mu = 0.0; beta = 0.5

# nearest neighbour coupling stencil
stencil = np.array([[0,1,0],[1,0,1],[0,1,0]])

#######################################################################

# split Metropolis-Hastings MCMC update
def step(i):
	global sigma
	force = J*convolve(sigma, stencil, mode='wrap') + mu*field
	trial = choice([-1.0,1.0], (n,n)); alpha = np.exp(beta*force*(trial-sigma))
	sigma = np.where((rand(n,n) > alpha) | ((parity+i)%2 == 0), sigma, trial)
	return sigma

#######################################################################

import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.lines import Line2D

fig = plt.figure(figsize=(48/5,27/5), frameon=False)
ax = fig.gca(); ax.set_aspect('equal')

# spin orientation map
crt = colors.LinearSegmentedColormap.from_list("CRT", ['black', 'greenyellow'])
spins = plt.imshow(sigma, origin='lower', vmin=-1.0, vmax=1.0, cmap=crt, interpolation='none')

# total magnetization meter
gauge = colors.LinearSegmentedColormap.from_list("magnet", ['blue', 'darkgreen', 'red'])
magnet = ax.add_line(Line2D([1.03*n,1.03*n], [n/2,n/2], color=gauge(0), linewidth=7, solid_capstyle='butt', clip_on=False, zorder=3))

# no ticks and tight framing
plt.tick_params(left=False, right=False, labelleft=False, labelbottom = False, bottom=False)
plt.tight_layout()

#######################################################################

# interactive controls
from matplotlib.widgets import Slider

# make room for widgets
fig.subplots_adjust(bottom=0.12)

# user interface elements
beta_slider = Slider(
    ax=fig.add_axes([0.1, 0.07, 0.8, 0.03]),
    valmin=0.0, valmax=1.0, valinit=beta,
    label='β'
)

mu_slider = Slider(
    ax=fig.add_axes([0.1, 0.02, 0.8, 0.03]),
    valmin=-1.0, valmax=1.0, valinit=mu,
    label='μ'
)

# update simulation parameters
def update(value):
	global beta, mu
	beta = beta_slider.val; mu = mu_slider.val

# register update handler
beta_slider.on_changed(update)
mu_slider.on_changed(update)

#######################################################################

import matplotlib.animation as animation

# called to advance animation to next frame
def animate(i):
	spins.set_data(step(i))
	M = (1.0 + np.sum(sigma)/n**2)/2.0
	magnet.set_data([1.03*n,1.03*n], [n/2,M*n])
	magnet.set_color(gauge(M))

animation = animation.FuncAnimation(fig, animate, frames=1000, interval=1000.0/60)
#animation.save('ising.mp4')

plt.show()
