#!/usr/bin/env python
# Floquet exponents and stability diagram (vectorized version)

#######################################################################

import numpy as np
from gl import gl12

#######################################################################

# plot and integrator resolution
n = 512; steps = 32

# period of the driver
t = np.pi; dt = t/steps

# parameter range to scan
a = np.linspace(-5.0,5.0,n)
q = np.linspace(0.0,10.0,n)

# batch parameters (passed as globals)
kappa = np.zeros(2*n)
twoq = 2.0*np.concatenate((q,q))

# vectorized dynamical system
def f(state):
	# principal solutions and driving oscillator
	x,v = state[:4*n].reshape(2,2*n); w,u = state[-2:]
	return np.concatenate((v,-(kappa-w*twoq)*x,np.array([u,-4.0*w])))

# compute Floquet exponents
def mu(state):
	trace = (state[0:n] + state[3*n:4*n])/2.0
	return np.arccosh(trace.astype(complex))/t

# compute a batch of Floquet exponents
def batch(a):
	# a is passed in a global variable 
	global kappa; kappa = np.full(2*n,a)
	# vectorized initial conditions
	state = np.concatenate((np.ones(n),np.zeros(2*n),np.ones(n),np.array([1.0,0.0])))
	for i in range(0,steps): state = gl12(f, state, dt)
	return np.real(mu(state))

#######################################################################

from time import perf_counter as now
from joblib import Parallel, delayed

t1 = now()
#floquet = np.array([batch(v) for v in a])
floquet = np.array(Parallel(n_jobs=8)(delayed(batch)(v) for v in a))
t2 = now(); print(t2-t1)

#######################################################################

'''
# fast decimator using CIC filter (for oversampled rendering)
# https://en.wikipedia.org/wiki/Cascaded_integrator–comb_filter
floquet = np.diff(np.cumsum(floquet, axis=0)[::4,:], axis=0)/4.0
floquet = np.diff(np.cumsum(floquet, axis=1)[:,::4], axis=1)/4.0
'''

#######################################################################

import matplotlib.pyplot as plt
import matplotlib.cm as cm

cmap = cm.afmhot; cmap.set_under('lightgray')

plt.imshow(floquet-(floquet==0.0), origin='lower', extent=[q[0],q[-1],a[0],a[-1]], vmin=0.0, cmap=cmap, norm='linear', aspect='equal', interpolation='none')
plt.colorbar()
plt.show()
