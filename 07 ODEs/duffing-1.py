#!/usr/bin/env python
# non-linear resonance in Duffing oscillator

#######################################################################

from math import *
import numpy as np
from gl import gl12

#######################################################################

# plot and integrator resolution
periods = 256; steps = 64

# parameters
alpha = 0.02
kappa = 1.0
beta = 0.1
force = 4.0
omega = 1.0

'''
alpha = 0.02
kappa = 1.0
beta = 5.0
force = 8.0
omega = 0.5
'''

#######################################################################

# period of the driver
t = 2.0*np.pi/omega; dt = t/steps

# dynamical system
def f(state):
	x,v,t = state
	return np.array([v,force*cos(t)-alpha*v-kappa*x-beta*x**3,omega])

# initial state
state = np.array([1.0,0.0,0.0])

#######################################################################

from time import perf_counter as now

# initialize Poincare map storage
history = np.zeros([periods,3])

# evolve the solution
t1 = now()
for k in range(periods):
	for i in range(steps):
		state = gl12(f, state, dt)
	history[k] = state
t2 = now(); print(t2-t1)

#######################################################################

import matplotlib.pyplot as plt

fig = plt.figure(); ax = plt.gca()

plt.plot(history[:,0], history[:,1], '.')

plt.xlim([-7,7])
plt.ylim([-7,7])
ax.set_aspect('equal')

plt.show()

