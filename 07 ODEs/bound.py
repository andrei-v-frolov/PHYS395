#!/usr/bin/env python
# animation of a bound orbit around a black hole

#######################################################################

from math import *
import numpy as np
from gl import gl10

#######################################################################

# minimal and maximal orbit radia
m = 0.5; rg = 2.0*m; a = 16.0; b = 64.0

# slowest timescale involved
T = 2.0*pi*sqrt(b**3/m)

# corresponding angular momentum
L2 = rg*a*a*b*b/((b-rg)*(a+b)*a - b*b*rg); L = sqrt(L2)

# equations of motion for a particle around BH
def f(state):
	r,pi,phi = state
	return [pi, (L2/r*(1.0-3.0*m/r)-m)/(r*r), L/r**2]

# energy (should be conserved for Hamiltonian EoM)
def E(state):
	r,pi,phi = state
	return pi*pi + (1.0-rg/r)*(1.0+(L/r)**2)

# initial conditions and energy
state = np.array([b,0.0,L/b**2]); E0 = E(state)

#######################################################################

import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.collections import LineCollection as lines

# plot the trajectory
fig = plt.figure(); ax = fig.gca()

# perihelion and aphelion
ax.add_patch(plt.Circle((0, 0), a, color='r', fill=False))
ax.add_patch(plt.Circle((0, 0), b, color='r', fill=False))

# black hole horizon
ax.add_patch(plt.Circle((0, 0), 1.0, color='black', fill=True))

# monitor energy violation
#deltaE = ax.text(0.0, 0.8*b, '', horizontalalignment='center', fontfamily='monospace')

# colormap interpolating color to full transparency
def fade(color, name=None):
	r,g,b,a = colors.to_rgba(color)
	return colors.LinearSegmentedColormap.from_list(name, [color, (r,g,b,0)])

# orbital trace fading away
trace = lines([], cmap=fade('skyblue'), linewidth=3)
tail = 1500; trace.set_array(np.linspace(0,1,tail))
ax.add_collection(trace)

# particle drawing (nothing fancy)
particle, = plt.plot([0], [0], "o", linewidth=7, ms=7)

ax.set_aspect('equal')
plt.xlim([-1.05*b,1.05*b])
plt.ylim([-1.05*b,1.05*b])

#######################################################################

import matplotlib.animation as animation

# initialize trace to starting position
history = np.empty([tail,2])
history[:,0].fill(state[0]*cos(state[2]))
history[:,1].fill(state[0]*sin(state[2]))

# called to advance animation to next frame
def animate(i):
	global state, history
	state = gl10(f, state, T/256)
	x = state[0]*cos(state[2])
	y = state[0]*sin(state[2])
	history = np.roll(history,1,axis=0); history[0] = [x,y]
	pts = history.reshape(-1, 1, 2)
	particle.set_data([x], [y])
	trace.set_paths(np.concatenate([pts[:-1], pts[1:]], axis=1))
	#deltaE.set_text('δE = %+0.1E' % (E(state)-E0))

animation = animation.FuncAnimation(fig, animate, frames=3000, interval=1000.0/60)
#animation.save('pendulum-2.mp4')
plt.show()
