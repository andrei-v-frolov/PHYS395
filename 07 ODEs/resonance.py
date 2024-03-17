#!/usr/bin/env python
# animated non-linear and parametric resonance demo

#######################################################################

from math import *
import numpy as np
from gl import gl10

#######################################################################

# simulation parameters
alpha = 0.1
omega = 1.0
kappa = 0.0
force = 0.0

# dynamical system to be integrated
def f(state):
	x,v,t = state
	return np.array([v,force*sin(t)-sin(x)/(1.0+kappa*sin(t))-alpha*v,omega])

# energy (should be conserved for Hamiltonian EoM)
def E(state):
	x,v,t = state
	return v*v/2.0 - cos(x)

# initial conditions and energy
state = np.array([0.0,0.0,0.0]); E0 = E(state)

#######################################################################

import matplotlib.pyplot as plt
from matplotlib.patches import Arc
from matplotlib.lines import Line2D

fig = plt.figure(); ax = fig.gca()

# circle traced by pendulum
ax.add_patch(plt.Circle((0, 0), 1.0, color='r', fill=False))

# arc representing driving force
arc = ax.add_patch(Arc((0, 0), 2.0, 2.0, color='orange', linewidth=7, zorder=3))

# line representing length modulation
line = ax.add_line(Line2D([0,0], [0,0], color='orange', linewidth=7, solid_capstyle='butt', zorder=3))

# total energy meter
energy = ax.add_line(Line2D([1,1], [-1,-1], color='green', linewidth=7, solid_capstyle='butt', clip_on=False, zorder=3))

# pendulum drawing (nothing fancy)
pendulum, = plt.plot([0,sin(state[0])], [0,-cos(state[0])], "o-", linewidth=7, ms=15)

ax.set_aspect('equal')
plt.xlim([-1.1,1.1])
plt.ylim([-1.1,1.1])

#######################################################################

# interactive controls
from matplotlib.widgets import Slider, RadioButtons

# make room for widgets
fig.subplots_adjust(bottom=0.22)

# user interface elements
freq_slider = Slider(
    ax=fig.add_axes([0.2, 0.12, 0.51, 0.03]),
    valmin=-2.0, valmax=2.0, valinit=0.0,
    label='log₂(ω/ω₀)'
)

value_slider = Slider(
    ax=fig.add_axes([0.2, 0.07, 0.51, 0.03]),
    valmin=0.0, valmax=1.0, valinit=0.0,
    label='amount'
)

alpha_slider = Slider(
    ax=fig.add_axes([0.2, 0.02, 0.51, 0.03]),
    valmin=0.0, valmax=0.3, valinit=alpha,
    label='losses'
)

kind = RadioButtons(
	fig.add_axes([0.8, 0.03, 0.18, 0.12]),
	('external','parametric'),
)

# update simulation parameters
def update(value):
	global alpha, omega, kappa, force
	alpha = alpha_slider.val
	omega = exp(log(2.0)*freq_slider.val)
	match kind.value_selected:
		case 'external':
			force = value_slider.val; kappa = 0.0
		case 'parametric':
			kappa = value_slider.val; force = 0.0
		case _: ...

# register update handler
freq_slider.on_changed(update)
value_slider.on_changed(update)
alpha_slider.on_changed(update)
kind.on_clicked(update)

#######################################################################

import matplotlib.animation as animation

# called to advance animation to next frame
def animate(i):
	global state
	state = gl10(f, state, 0.03)
	x,v,t = state; s = sin(x); c = cos(x)
	pendulum.set_data([0,s], [0,-c])
	energy.set_data([1.1,1.1], [-1,E(state)])
	line.set_data([s,(1.0+kappa*sin(t))*s], [-c,-(1.0+kappa*sin(t))*c])
	a = (180.0/pi)*x - 90.0; b = a + (90.0/pi)*force*sin(t)
	arc.theta1 = min(a,b); arc.theta2 = max(a,b)

animation = animation.FuncAnimation(fig, animate, frames=1000, interval=1000.0/60)
plt.show()
