#!/usr/bin/env python
# time it takes for double pendulum to flip (fancy version)
# run as: python dpflip.py [{save|resume|plot} file.npz]

#######################################################################

import numpy as np
from gl import gl12

#######################################################################

# plot resolution and integration length
n = 64; dt = np.pi/32; steps = 320; epochs = 10

# optimal vectorized chunk and # of threads
chunk = 1024; jobs = 8

# parameter range to scan
a = np.pi * (n-1)/n
x = np.linspace(-a,a,n)
y = np.linspace(-a,a,n)
X,Y = np.meshgrid(x,y)

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

# all initial momenta are zero
P = np.zeros(n*n)

# initial conditions to scan (flattened)
state = np.array((X.flat,Y.flat,P,P)); E0 = E(state)

# ICs with enough energy to flip (flattened)
mask = (3.0*np.cos(X) + np.cos(Y) < 2.0).reshape(-1)

# time when flipped (negative if not flipped yet)
T = -np.ones(n*n); epoch = 0.0

#######################################################################

# access command line arguments
from sys import argv

# default checkpoint file name
savefile = f'dpflip-{n}.npz'

# save current evolution state to a file
def checkpoint(file):
	np.savez_compressed(file, time=T, state=state, epoch=epoch)

# load previously saved evolution state
def resume(file):
	global T, state, mask, epoch; data = np.load(file)
	T = data['time']; epoch = data['epoch']
	state = data['state']; mask &= (T < 0.0)

# parse command line arguments
if len(argv) > 1:
	match argv[1:]:
		case ('save', file):
			savefile = file
		case ('resume', file):
			savefile = file
			resume(file)
		case ('plot', file):
			epochs = -1
			resume(file)
		case _:
			...

#######################################################################

from time import perf_counter as now
from joblib import Parallel, delayed

# evolve a batch of (vectorized) pendulum instances
def batch(state):
	# time of the first flip (negative if not flipped yet)
	t = -np.ones(state.shape[-1])
	
	# evolve and record the time of first flip
	for i in range(steps):
		state = gl12(f, state.reshape(-1), dt).reshape(4,-1)
		t = np.where((t < 0.0) & flipped(state), epoch + (i+1)*dt, t)

	# return stacked results (for easier parallel merge)
	return np.vstack((state,t))

# dispatch workload for parallel evolution (in batches)
def evolve(state):
	batches = state.shape[-1]//chunk + 1
	workload = np.array_split(state, batches, axis=-1)
	results = Parallel(n_jobs=jobs)(delayed(batch)(s) for s in workload)
	return np.hstack(results)

# evolve non-flipped pendulae, check-pointing evolution
for i in range(epochs):
	print(f'epoch {i}: evolving {mask.sum()} pixels', end='')
	t1 = now()
	data = evolve(state[:,mask])
	state[:,mask] = data[:4]; T[mask] = data[-1]
	t2 = now(); print(f', {t2-t1:.3f}s')
	epoch += steps*dt; mask &= (T < 0.0)
	checkpoint(savefile)

# plot conserved energy violation instead
#T = E(state)-E0

# reshape the result to original grid size
T = T.reshape(n,n)

#######################################################################
'''
from scipy.ndimage import median_filter as median

# downsample using median filter
T = median(T, 3)[1::3,1::3]
'''
#######################################################################

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm

cmap = cm.turbo; cmap.set_bad('lightgray')
slog = colors.SymLogNorm(linthresh=1.0e-8, vmin=-1.0e-4, vmax=1.0e-4)

plt.imshow(T, origin='lower', extent=[x[0],x[-1],y[0],y[-1]], cmap=cmap, vmin=dt, norm='log', aspect='equal', interpolation='none')
#plt.imshow(T, origin='lower', extent=[x[0],x[-1],y[0],y[-1]], cmap='RdBu_r', norm=slog, aspect='equal', interpolation='none')
plt.colorbar()

#plt.contourf(X, Y, 3.0*np.cos(X) + np.cos(Y), [2.0,5.0], colors='lightgray')

plt.show()
