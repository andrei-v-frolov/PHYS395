#!/usr/bin/env python
# time it takes for double pendulum to flip

#######################################################################

import numpy as np
from gl import gl12

#######################################################################

# plot resolution and integration length
n = 64; dt = np.pi/32; steps = 320; jobs = 8

# parameter range to scan
a = np.pi * (1.0 - 1.0/n)
x = np.linspace(-a,a,n); y = x

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
def batch(value, x1=x):
	# vectorized initial conditions
	t = np.zeros_like(x1)
	x2 = np.broadcast_to(value,x1.shape)
	state = np.concatenate((x1,x2,t,t)); t = -1.0

	# record the time of first flip
	for i in range(0,steps):
		state = gl12(f, state, dt)
		t = np.where((t < 0.0) & flipped(state), (i+1)*dt, t)
	return t

#######################################################################

from time import perf_counter as now
from joblib import Parallel, delayed

t1 = now(); print("Scanning %ix%i IC grid:" % (n,n))
T = np.array(Parallel(n_jobs=jobs)(delayed(batch)(theta) for theta in y))
t2 = now(); print(t2-t1)

# checkpoint storing data
np.save('dpflip.npy', T)

#######################################################################

# rescan ICs with enough energy to flip
X,Y = np.meshgrid(x,y); steps *= 10
mask = (T < 0.0) & (3.0*np.cos(X) + np.cos(Y) < 2.0)

# split the workload into batches
x1 = np.array_split(np.extract(mask,X),jobs)
x2 = np.array_split(np.extract(mask,Y),jobs)

t1 = now(); print("Rescanning %i pixels:" % np.sum(mask))
details = np.concatenate(Parallel(n_jobs=jobs)(delayed(batch)(i,j) for i,j in zip(x2,x1)))
t2 = now(); print(t2-t1)

# put them back in place
np.place(T, mask, details)

# checkpoint storing data
np.save('dpflip.npy', T)

#######################################################################
'''
from scipy.ndimage import median_filter as median

# load checkpoint data
#T = np.load('dpflip.npy')

# downsample using median filter
T = median(T, 3)[1::3,1::3]
'''
#######################################################################

import matplotlib.pyplot as plt
import matplotlib.cm as cm

cmap = cm.turbo; cmap.set_bad('lightgray')

plt.imshow(T, origin='lower', extent=[x[0],x[-1],y[0],y[-1]], cmap=cmap, vmin=dt, norm='log', aspect='equal')
plt.colorbar()
plt.show()
