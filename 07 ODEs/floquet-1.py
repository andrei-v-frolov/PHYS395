#!/usr/bin/env python
# parametric resonance and Floquet exponents

#######################################################################

import numpy as np
from gl import gl12

#######################################################################

# integration length and resolution
periods = 16; steps = 64; n = periods*steps

# period of the driver
t = np.pi; dt = t/steps

# parameter values
a = 1.0; q = 0.10

# dynamical system
def f(state):
	x1,x2,v1,v2,w,u = state; kappa = a - 2.0*q*w
	return np.array([v1,v2,-kappa*x1,-kappa*x2,u,-4.0*w])

# compute Floquet exponent
def mu(state):
	x1,x2,v1,v2,w,u = state; trace = (x1+v2)/2.0
	return np.arccosh(trace.astype(complex))/(n*dt)

# initial conditions
state = np.array([1.0,0.0,0.0,1.0,1.0,0.0])

#######################################################################

from time import perf_counter as now

history = np.zeros([n,6])

t1 = now()
# evolve the dynamical system
for i in range(n):
	state = gl12(f, state, dt)
	history[i] = state
t2 = now(); print(t2-t1)

# Floquet exponent (evaluate over integer number of periods!)
print(mu(state))

#######################################################################

import matplotlib.pyplot as plt

p = np.arange(n)/steps
plt.plot(p, history[:,4], '-')
plt.plot(p, history[:,1], 'r-', linewidth=3.0)
plt.xlim([0,n/steps])

plt.show()
