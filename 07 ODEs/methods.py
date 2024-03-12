#!/usr/bin/env python
# demo of ODE integration methods

#######################################################################

import numpy as np
from scipy.optimize import root

#######################################################################

# dimension of the state vector
dim = 2

# number and length of steps to take
n = 1024; dt = 10.0*np.pi/n

# dynamical system to be integrated
def f(state):
	x,v = state
	return np.array([v,-x**3])

# energy (should be conserved for Hamiltonian EoM)
def E(state):
	x,v = state
	return v*v/2.0 + x**4/4.0

# initial conditions and energy
state = np.array([1.0,0.0]); E0 = E(state)

#######################################################################
# Euler and Runge-Kutta methods
#######################################################################

# forward Euler step (1st order)
def euler(state):
	return state + f(state)*dt

# backward Euler step (1st order)
def ieuler(state):
	return root(lambda x: x - state - f(x)*dt, state).x

# explicit midpoint (aka RK2)
def midpoint(state):
	return state + f(state + f(state)*dt/2.0)*dt

# implicit midpoint (2nd order)
def imidpoint(state):
	return root(lambda x: x - state - f((x+state)/2.0)*dt, state).x

# to be continued...
...

#######################################################################
# operator splitting methods for separable Hamiltonian (hard-coded EoM)
#######################################################################

# 1st order Hamiltonian split
def si1(state):
	x,v = state
	x += v*dt
	v += -x**3*dt
	return np.array([x,v])

# 2nd order Hamiltonian split
def si2(state):
	x,v = state
	x += v*dt/2.0
	v += -x**3*dt
	x += v*dt/2.0
	return np.array([x,v])

# to be continued...
...

#######################################################################
# Gauss-Legendre methods; symplectic with arbitrary Hamiltonian, A-stable
#######################################################################

...

#######################################################################

# evolution history and violation of energy conservation
t = 0.0; history = np.zeros([dim+2,n])

# evolve dynamical system with specified method
for i in range(0,n):
	t += dt; state = midpoint(state)
	history[:,i] = [*state, t, E(state)-E0]

#######################################################################

import matplotlib.pyplot as plt

#plt.plot(history[dim], history[0])
plt.plot(history[dim], history[dim+1], "r-")
plt.show()
