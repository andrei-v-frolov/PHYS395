#!/usr/bin/env python
# animated physical pendulum evolution

#######################################################################

from math import *
import numpy as np
from gl import gl10

#######################################################################

# dynamical system to be integrated
def f(state):
	x,v = state
	return np.array([v,-sin(x)])

# energy (should be conserved for Hamiltonian EoM)
def E(state):
	x,v = state
	return v*v/2.0 - cos(x)

# initial conditions and energy
state = np.array([1.0,0.0]); E0 = E(state)

#######################################################################

'''
from scipy.integrate import solve_ivp

# stopping condition (half period completed)
def stop(t,state):
	x,v = state
	return v

# stopping condition attributes
stop.terminal = True
stop.direction = 1.0

# edit f(y) -> f(t,y) to use this solver
soln = solve_ivp(f, [0.0,6.0], state, events=stop, rtol=1e-16, atol=1e-16)

period = 2.0*soln.t[-1]; print(period); exit()
'''

#######################################################################

import matplotlib.pyplot as plt

fig = plt.figure(); ax = fig.gca()

# circle traced by pendulum
ax.add_patch(plt.Circle((0, 0), 1.0, color='r', fill=False))

# monitor energy violation
#deltaE = ax.text(0.0, 0.8, '', horizontalalignment='center', fontfamily='monospace')

# pendulum drawing (nothing fancy)
pendulum, = plt.plot([0,sin(state[0])], [0,-cos(state[0])], "o-", linewidth=7, ms=15)

ax.set_aspect('equal')
plt.xlim([-1.1,1.1])
plt.ylim([-1.1,1.1])

#######################################################################

import matplotlib.animation as animation

# called to advance animation to next frame
def animate(i):
	global state
	state = gl10(f, state, 6.699975664370446/200); x = state[0]
	pendulum.set_data([0,sin(x)], [0,-cos(x)])
	#deltaE.set_text('δE = %+0.1E' % (E(state)-E0))

animation = animation.FuncAnimation(fig, animate, frames=1000, interval=1000.0/60)
#animation.save('pendulum.mp4')
plt.show()
