#!/usr/bin/env python
# demo of ODE integration methods

#######################################################################

from math import *
import numpy as np
from scipy.optimize import root

#######################################################################
# anharmonic oscillator (two-dimensional dynamic system example)
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

'''
#######################################################################
# Kepler problem (four-dimensional dynamic system example)
#######################################################################

# dimension of the state vector
dim = 4

# number and length of steps to take
n = 2048; dt = 2.0*np.pi/16

# dynamical system to be integrated
def f(state):
	x,y,vx,vy = state; r = sqrt(x*x+y*y)
	return np.array([vx,vy,-x/r**3,-y/r**3])

# energy (should be conserved for Hamiltonian EoM)
def E(state):
	x,y,vx,vy = state; r = sqrt(x*x+y*y)
	return (vx*vx+vy*vy)/2.0 - 1.0/r

# initial conditions and energy
state = np.array([1.0,0.0,0.0,-1.0]); E0 = E(state)
'''

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

# implicit midpoint (aka GL2)
def imidpoint(state):
	return root(lambda x: x - state - f((x+state)/2.0)*dt, state).x

# 4-th order Runge-Kutta method
def rk4(state):
	k1 = f(state)
	k2 = f(state+k1*dt/2.0)
	k3 = f(state+k2*dt/2.0)
	k4 = f(state+k3*dt)
	return state + (k1+2.0*k2+2.0*k3+k4)*dt/6.0

#######################################################################
# operator splitting methods for separable Hamiltonian (hard-coded EoM)
#######################################################################

# 1st order Hamiltonian split (oscillator)
def si1(state):
	x,v = state
	x += v*dt
	v += -x**3*dt
	return np.array([x,v])

# 2nd order Hamiltonian split (oscillator)
def si2(state):
	x,v = state
	x += v*dt/2.0
	v += -x**3*dt
	x += v*dt/2.0
	return np.array([x,v])

'''
# 2nd order Hamiltonian split (Kepler problem)
def si2(state):
	x,y,vx,vy = state
	x += vx*dt/2.0
	y += vy*dt/2.0
	r = sqrt(x*x+y*y)
	vx -= x/r**3*dt
	vy -= y/r**3*dt
	y += vy*dt/2.0
	x += vx*dt/2.0
	return np.array([x,y,vx,vy])
'''

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
	t += dt; state = rk4(state)
	history[:,i] = [*state, t, E(state)-E0]

#######################################################################

import matplotlib.pyplot as plt

#plt.plot(history[dim], history[0])
plt.plot(history[dim], history[dim+1], "r-")

'''
#plt.plot(history[0], history[1])
plt.gca().set_aspect('equal')
plt.xlim([-1.5,1.5])
plt.ylim([-1.5,1.5])
'''

plt.show()
