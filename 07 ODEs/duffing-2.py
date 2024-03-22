#!/usr/bin/env python
# non-linear resonance in Duffing oscillator (vectorized version)

#######################################################################

from math import *
import numpy as np
from gl import gl12

#######################################################################

# plot and integrator resolution
n = 1000; periods = 256; steps = 32

# oscillator parameters
alpha = 0.2
kappa = -1.0
beta  = 1.0
force = 0.3
omega = 1.0

'''
alpha = 0.02
kappa = 1.0
beta  = 5.0
force = 8.0
omega = 0.5
'''

#######################################################################

# period of the driver
t = 2.0*np.pi/omega; dt = t/steps

# vectorized dynamical system
def f(state):
	x,v,t = state.reshape(3,len(state)//3)
	return np.concatenate((v,force*np.cos(t)-alpha*v-(kappa+beta*x*x)*x,omegas))

# initial state placeholder
state = np.zeros(3*n)

# seed the initial state evolving over the driver period
omegas = np.array([omega]); ic = np.array([2.0,0.0,0.0])

# propagate the seed IC to batch initial state vector
for i in range(n):
	ic = gl12(f, ic, t/n)
	state[i::n] = ic

#######################################################################

from time import perf_counter as now

# initialize Poincare map storage
omegas = np.full(n,omega)
history = np.zeros([periods,3*n])

# evolve the initial conditions batch
t1 = now()
for k in range(periods):
	for i in range(steps):
		state = gl12(f, state, dt)
	history[k] = state
t2 = now(); print(t2-t1)

#######################################################################

import matplotlib.pyplot as plt

fig = plt.figure(); ax = fig.gca()

poincare = plt.scatter(history[:,0], history[:,n], marker=".", s=1.0)
plt.xlim([-2,2])
plt.ylim([-1,1])

#######################################################################

import matplotlib.animation as animation

# called to advance animation to next frame
def animate(i):
	poincare.set_offsets(history[:,i::n])

animation = animation.FuncAnimation(fig, animate, frames=n, interval=1000.0/60)
#animation.save('poincare.mp4')
plt.show()
