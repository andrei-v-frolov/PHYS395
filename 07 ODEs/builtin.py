#!/usr/bin/env python
# demo of built-in SciPy integrators

#######################################################################

import numpy as np
from scipy.integrate import solve_ivp

#######################################################################

# number of steps to take
n = 1024; tmax = 10.0*np.pi

# dynamical system to be integrated
def f(t,state):
	x,v = state; return np.array([v,-x**3])

# energy (should be conserved for Hamiltonian EoM)
def E(state):
	x,v = state; return v*v/2.0 + x**4/4.0

# initial conditions and energy
state = np.array([1.0,0.0]); E0 = E(state)

#######################################################################

# time range to evolve for
t = np.linspace(0.0, tmax, n)

'''
# stopping condition (add events=stop to driver)
def stop(t,state):
	x,v = state; return v

# stopping condition attributes (e.g. terminate)
stop.terminal = True
stop.direction = 1.0
'''

soln = solve_ivp(f, [t[0],t[-1]], state, method='RK45', t_eval=t)

#######################################################################

import matplotlib.pyplot as plt

#plt.plot(soln.t, soln.y[0])
plt.plot(soln.t, np.apply_along_axis(E,0,soln.y)-E0, "r-")
plt.show()
