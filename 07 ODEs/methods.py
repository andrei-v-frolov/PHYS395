#!/usr/bin/env python
# demo of ODE integration methods

#######################################################################

from math import *
import numpy as np
from scipy.optimize import root

#######################################################################
# anharmonic oscillator (two-dimensional dynamic system example)
#######################################################################

# number and length of steps to take
n = 1024; dt = 10.0*np.pi/n

# dynamical system to be integrated
def f(state):
	x,v = state; return np.array([v,-x**3])

# energy (should be conserved for Hamiltonian EoM)
def E(state):
	x,v = state; return v*v/2.0 + x**4/4.0

# initial conditions and energy
state = np.array([1.0,0.0]); E0 = E(state)

#######################################################################
# Kepler problem (six-dimensional dynamic system example)
#######################################################################
'''
# number and length of steps to take
n = 2048; dt = 2.0*np.pi/16

# dynamical system to be integrated
def f(state):
	x,v = state.reshape([2,3]); r = sqrt(sum(x*x))
	return np.concatenate([v,-x/r**3])

# energy (should be conserved for Hamiltonian EoM)
def E(state):
	x,v = state.reshape([2,3]); r = sqrt(sum(x*x))
	return sum(v*v)/2.0 - 1.0/r

# initial conditions and energy
state = np.array([1.0,0.0,0.0, 0.0,1.0,0.0]); E0 = E(state)
'''

#######################################################################
# Euler and Runge-Kutta methods
#######################################################################

# forward Euler step (1st order)
def euler(state, dt):
	return state + f(state)*dt

# backward Euler step (1st order)
def ieuler(state, dt):
	return root(lambda x: state + f(x)*dt - x, state).x

# explicit midpoint (aka RK2)
def midpoint(state, dt):
	return state + f(state + f(state)*dt/2.0)*dt

# implicit midpoint (aka GL2)
def imidpoint(state, dt):
	return root(lambda x: state + f((state+x)/2.0)*dt - x, state).x

# 4-th order Runge-Kutta method
def rk4(state, dt):
	k1 = f(state)
	k2 = f(state + k1*dt/2.0)
	k3 = f(state + k2*dt/2.0)
	k4 = f(state + k3*dt)
	return state + (k1 + 2.0*k2 + 2.0*k3 + k4)*dt/6.0

# n-th order Richardson extrapolation (even n only)
def re(n, state, dt):
	match n:
		case 1: return euler(state, dt)
		case 2: return midpoint(state, dt)
		case 4: return rk4(state, dt)
		case k:
			w = 2**(k-2)
			y1 = re(k-2, state, dt)
			y2 = re(k-2, state, dt/2.0)
			y3 = re(k-2, y2, dt/2.0)
			return (w*y3 - y1)/(w - 1.0)

#######################################################################
# operator splitting methods for separable Hamiltonian (hard-coded EoM)
#######################################################################

# 1st order Hamiltonian split (oscillator)
def si1(state, dt):
	x,v = state
	x +=  v*dt
	v += -x**3*dt
	return np.array([x,v])

# 2nd order Hamiltonian split (oscillator)
def si2(state, dt):
	x,v = state
	x +=  v*dt/2.0
	v += -x**3*dt
	x +=  v*dt/2.0
	return np.array([x,v])

'''
# 2nd order Hamiltonian split (Kepler problem)
def si2(state, dt):
	x,v = state.reshape([2,3])
	x +=  v*dt/2.0
	r = sqrt(sum(x*x))
	v += -x/r**3*dt
	x +=  v*dt/2.0
	return np.concatenate([x,v])
'''

# Yoshida 6-th order scheme timesteps
W6 = np.array([
	 1.31518632068391121888424972823886251E0,
	-1.17767998417887100694641568096431573E0,
	 0.235573213359358133684793182978534602E0,
	 0.784513610477557263819497633866349876E0
])

# 6th order symplectic integrator
def si6(state, dt):
	for i in range(-3,4):
		state = si2(state, W6[abs(i)]*dt)
	return state

# n-th order symplectic integrator (even n only)
def si(n, state, dt):
	match n:
		case 1: return si1(state, dt)
		case 2: return si2(state, dt)
		case 6: return si6(state, dt)
		case k:
			gamma = 1.0/(k-1.0)
			alpha = 1.0/(2.0 - 2.0**gamma)
			state = si(k-2, state, alpha*dt)
			state = si(k-2, state, (1.0-2.0*alpha)*dt)
			state = si(k-2, state, alpha*dt)
			return state

#######################################################################
# Gauss-Legendre methods; symplectic with arbitrary Hamiltonian, A-stable
#######################################################################

from gl import gl4, gl6, gl8, gl10, gl12

#######################################################################

# evolution history and violation of energy conservation
t = np.arange(1,n+1)*dt; history = np.zeros([n,len(state)+1])

# evolve dynamical system with specified method
for i in range(n):
	state = rk4(state, dt)
	#state = gl4(f, state, dt)
	history[i] = [*state, E(state)-E0]

#######################################################################

import matplotlib.pyplot as plt

#plt.plot(t, history)
plt.plot(t, history[:,-1], 'r-')

'''
plt.plot(history[:,0], history[:,1], '-')
plt.gca().set_aspect('equal')
plt.xlim([-1.5,1.5])
plt.ylim([-1.5,1.5])
'''

plt.show()

