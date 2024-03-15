#!/usr/bin/env python
# plot exact and integrated phase space trajectories

#######################################################################

import numpy as np

#######################################################################

# physical pendulum EoM
def si2(state,dt):
	x,v = state
	x += v*dt/2.0
	v -= np.sin(x)*dt
	x += v*dt/2.0
	return [x,v]

# energy of physical pendulum
def H(x,v):
	return v*v/2.0 - np.cos(x)

# initial points to sample trajectories
vs = np.linspace(0.0,4.0,9)
levels = [H(0.0,v) for v in vs]

#######################################################################

# sample orbit of symplectic integrator
def orbit(x0,v0,dt=2.0*np.pi/16,n=1024):
	y = np.zeros([n,2]); y[0] = [x0,v0]
	for i in range(1,n): y[i] = si2(y[i-1],dt)
	idx = np.arctan2(y[:,1],y[:,0]).argsort()
	return y[idx,0], y[idx,1]

#######################################################################

import matplotlib.pyplot as plt

# phase space region to plot
n = 1024; xmax = 10.0; ymax = 4.5

# sample energy functional on a grid
xgrid = np.linspace(-xmax, xmax, n)
ygrid = np.linspace(-ymax, ymax, n)
X,Y = np.meshgrid(xgrid, ygrid)
E = H(X,Y)

fig = plt.figure()

# plot algebraic phase space trajectories
plt.contour(X,Y,E, levels=levels, colors='tab:blue')

# plot integrated phase space trajectories
for v in vs[:-4]:
	q,p = orbit(0.0,v)
	plt.fill(q,p,edgecolor='tab:red',linewidth=3,fill=False)

# set square aspect ratio
plt.gca().set_aspect('equal')
plt.xlim([-xmax,xmax])
plt.ylim([-ymax,ymax])

plt.show()
