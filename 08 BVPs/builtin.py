#!/usr/bin/env python
# solve soliton BVP using SciPy built-in solver

#######################################################################

import numpy as np
from scipy.integrate import solve_bvp

#######################################################################

# field theory Lagrangian parameters
alpha = 1.0; nu = 1.0

# dynamical system corresponding to EoM
def f(t, state):
	x,v = state
	return [v,alpha*(x*x-nu*nu)*x]

# boundary conditions
def bc(a,b):
	return [a[0],b[0]-nu]

# built-in solver
x = np.linspace(0.0, 100.0, 1024)
y = [x/(1.0+x*x)**0.5, (1.0+x*x)**(-1.5)]
soln = solve_bvp(f, bc, x, y, max_nodes=16*1024)

#######################################################################

import matplotlib.pyplot as plt

# plot sampled solution
plt.plot(soln.x, soln.y[0])
#plt.plot(soln.x, soln.y[0]-np.tanh(soln.x/np.sqrt(2.0)), 'r')
#plt.xlim([0.0,10.0])
#plt.ylim([-1e-7,1e-7])

plt.show()
