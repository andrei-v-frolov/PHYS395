#!/usr/bin/env python
# time it takes for double pendulum to flip

#######################################################################

import numpy as np
from gl import gl12

#######################################################################

# plot resolution and integration length
n = 64; dt = np.pi/32; steps = 320; jobs = 8

# parameter range to scan
x = np.linspace(-3.14,3.14,n); y = x

#######################################################################

# compute angular velocities
def v(state):
	x1,x2,p1,p2 = state.reshape(4,-1)
	sigma = np.cos(x1-x2)
	kappa = (16.0/9.0) - sigma*sigma
	a = (2.0/3.0)*p1 - sigma*p2
	b = (8.0/3.0)*p2 - sigma*p1
	return (a/kappa,b/kappa)

# dynamical system to be integrated
def f(state):
	x1,x2,p1,p2 = state.reshape(4,-1)
	v1,v2 = v(state); w = v1*v2 * np.sin(x1-x2)
	return np.concatenate((v1,v2,-3.0*np.sin(x1)-w,-np.sin(x2)+w))

# energy (should be conserved for Hamiltonian EoM)
def E(state):
	x1,x2,p1,p2 = state.reshape(4,-1); v1,v2 = v(state)
	K = (4.0*v1*v1+v2*v2)/3.0 + np.cos(x1-x2)*v1*v2
	return K-(3.0*np.cos(x1) + np.cos(x2))

# check if either of the arms flipped
def flipped(state):
	x1,x2,p1,p2 = state.reshape(4,-1)
	return (np.abs(x1) > np.pi) | (np.abs(x2) > np.pi)

#######################################################################

# compute an IC scanline in a batch
def batch(theta):
	# vectorized initial conditions
	state = np.concatenate((x,np.full(n,theta),np.zeros(2*n))); T = -np.ones(n)
	for i in range(0,steps):
		state = gl12(f, state, dt)
		T = np.where((T < 0.0) & flipped(state), (i+1)*dt, T)
	return T

#######################################################################

from time import perf_counter as now
from joblib import Parallel, delayed

t1 = now(); print("Scanning IC: ", end='')
T = np.array(Parallel(n_jobs=jobs)(delayed(batch)(theta) for theta in y))
t2 = now(); print(t2-t1)

#######################################################################

import matplotlib.pyplot as plt
import matplotlib.cm as cm

cmap = cm.turbo; cmap.set_bad('lightgray')

plt.imshow(T, origin='lower', extent=[x[0],x[-1],y[0],y[-1]], cmap=cmap, norm='log', aspect='equal', interpolation='none')
plt.colorbar()
plt.show()
