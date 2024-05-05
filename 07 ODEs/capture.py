#!/usr/bin/env python
# interactive demo of trajectories around a black hole

#######################################################################

from math import *
import numpy as np
from scipy.integrate import solve_ivp

#######################################################################

# momentum and impact parameter
pi = 1.0; b = 3.334; L = pi*b

# equations of motion for a particle around BH
def f(t,state):
	rho,pi = state
	return [pi,0.5-L*L*rho*(1.0-1.5*rho)]

# energy (should be conserved for Hamiltonian EoM)
def E(state):
	rho,pi = state
	return pi*pi + (1.0-rho)*(1.0+(L*rho)**2)

# stopping condition (paticle escapes or is captured)
def stop(t,state):
	rho,pi = state
	return rho*(1.0-rho)

# stopping condition attributes (terminate on escape)
stop.terminal = True
stop.direction = -1.0

#######################################################################

# integrate particle trajectory
def trajectory(pi, b, n=256):
	global L; L = pi*b; t = np.linspace(0.0, 10.0/(pi+1.0), n)
	soln = solve_ivp(f, [t[0],t[-1]], [0,pi], t_eval=t, events=stop, method='Radau')
	phi = L*soln.t[1:]; r = 1.0/soln.y[0,1:]
	return (-r*np.cos(phi),r*np.sin(phi))

x,y = trajectory(pi,b)

#######################################################################

import matplotlib.pyplot as plt

# plot the trajectory
fig = plt.figure(); ax = fig.gca()
orbit, = plt.plot(x, y, linewidth=3, zorder=0)

# black hole horizon
ax.add_patch(plt.Circle((0, 0), 1.0, color='black', fill=True))

# set maximal bounds
bmax = 10

ax.set_aspect('equal')
plt.xlim([-1.3*bmax,1.3*bmax])
plt.ylim([-bmax,bmax])

#######################################################################

# interactive controls
from matplotlib.widgets import Slider

# make room for widgets
fig.subplots_adjust(left=0.2, bottom=0.10)

# user interface elements
b_slider = Slider(
    ax=fig.add_axes([0.07, 0.13, 0.0225, 0.72]),
    valmin=-bmax, valmax=bmax, valinit=0,
    orientation="vertical",
    label="b"
)

pi_slider = Slider(
    ax=fig.add_axes([0.2, 0.02, 0.7, 0.03]),
    valmin=-5.0, valmax=5.0, valinit=0.0,
    label='log₂π'
)

# function to be called anytime a slider value changes
def update(value):
	x,y = trajectory(2**pi_slider.val, b_slider.val)
	orbit.set_data(x,y); #fig.canvas.draw_idle()

# register update handler
b_slider.set_val(b)
b_slider.on_changed(update)
pi_slider.on_changed(update)

plt.show()
