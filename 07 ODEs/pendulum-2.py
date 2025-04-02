#!/usr/bin/env python
# animated double pendulum evolution

#######################################################################

from math import *
import numpy as np
from gl import gl10

#######################################################################

# compute angular velocities
def v(state):
	x1,x2,p1,p2 = state
	sigma = cos(x1-x2)
	kappa = (16.0/9.0) - sigma*sigma
	a = (2.0/3.0)*p1 - sigma*p2
	b = (8.0/3.0)*p2 - sigma*p1
	return (a/kappa,b/kappa)

# dynamical system to be integrated
def f(state):
	x1,x2,p1,p2 = state; v1,v2 = v(state); w = v1*v2 * sin(x1-x2)
	return np.array([v1,v2,-3.0*sin(x1)-w,-sin(x2)+w])

# energy (should be conserved for Hamiltonian EoM)
def E(state):
	x1,x2,p1,p2 = state; v1,v2 = v(state)
	K = (4.0*v1*v1+v2*v2)/3.0 + cos(x1-x2)*v1*v2
	return K-(3.0*cos(x1) + cos(x2))

# initial conditions and energy
state = np.array([2*pi/3.0,2*pi/3.0,0.0,0.0]); E0 = E(state)

#######################################################################

import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.collections import LineCollection as lines

fig = plt.figure(); ax = fig.gca()

# circle traced by pendulum
ax.add_patch(plt.Circle((0, 0), 1.0, color='r', fill=False))
ax.add_patch(plt.Circle((0, 0), 2.0, color='r', fill=False))

# monitor energy violation
#deltaE = ax.text(0.0, 1.6, '', horizontalalignment='center', fontfamily='monospace')

# colormap interpolating color to full transparency
def fade(color, name=None):
	r,g,b,a = colors.to_rgba(color)
	return colors.LinearSegmentedColormap.from_list(name, [color, (r,g,b,0)])

# pendulum trace fading away
trace = lines([], cmap=fade('skyblue'), linewidth=3)
tail = 500; trace.set_array(np.linspace(0,1,tail))
ax.add_collection(trace)

# pendulum drawing (nothing fancy)
pendulum, = plt.plot([0,0,0], [0,0,0], "o-", linewidth=7, ms=15)

ax.set_aspect('equal')
plt.xlim([-2.1,2.1])
plt.ylim([-2.1,2.1])

#######################################################################

import matplotlib.animation as animation

# initialize trace to starting position
history = np.empty([tail,2])
history[:,0].fill( sin(state[0]) + sin(state[1]))
history[:,1].fill(-cos(state[0]) - cos(state[1]))

# called to advance animation to next frame
def animate(i):
	global state, history
	state = gl10(f, state, pi/100)
	x1 = sin(state[0]); x2 = x1 + sin(state[1])
	y1 = cos(state[0]); y2 = y1 + cos(state[1])
	history = np.roll(history,1,axis=0); history[0] = [x2,-y2]
	pts = history.reshape(-1, 1, 2)
	pendulum.set_data([0,x1,x2], [0,-y1,-y2])
	trace.set_paths(np.concatenate([pts[:-1], pts[1:]], axis=1))
	#deltaE.set_text('δE = %+0.1E' % (E(state)-E0))

animation = animation.FuncAnimation(fig, animate, frames=1000, interval=1000.0/60)
#animation.save('pendulum-2.mp4')
plt.show()
