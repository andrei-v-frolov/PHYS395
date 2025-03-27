#!/usr/bin/env python
# solve soliton BVP using shooting method

#######################################################################

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import bisect

#######################################################################

# field theory Lagrangian parameters
alpha = 1.0; nu = 1.0

# dynamical system corresponding to EoM
def f(t, state):
	x,v = state
	return [v,alpha*(x*x-nu*nu)*x]

# map IC at origin to BC violation far away
def bc(v0, infinity=10.0):
	soln = solve_ivp(f, [0,infinity], [0,v0])
	return soln.y[0,-1] - nu

# find IC satisfying desired BC (use robust solver!)
v0 = bisect(bc, 0.5, 1.0)

#######################################################################

import matplotlib.pyplot as plt

# sample found solution
t = np.linspace(0.0, 10.0, 1024)
soln = solve_ivp(f, [0,t[-1]], [0,v0], t_eval=t)

# plot sampled solution
#plt.plot(soln.t, soln.y[0])
plt.plot(soln.t, soln.y[0]-np.tanh(soln.t/np.sqrt(2.0)), 'r')

plt.show()
