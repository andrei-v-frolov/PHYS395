#!/usr/bin/env python
# cylindrical waveguide TM modes (spectral solver for radial eigenfunctions)

###############################################################################

import numpy as np
from scipy.linalg import eig, norm

###############################################################################

# grid resolution and waveguide radius
n = 256; r = 1.0

# wavenumbers of mode to plot
l = 4; m = 2

# initialize 1D rho & phi grids
twopi = 2.0*np.pi; halfpi = np.pi/2.0; dtheta = halfpi/n

theta = np.linspace(halfpi-dtheta/2.0, dtheta/2.0, n)
rho = r*np.cos(theta); phi = np.linspace(0.0, twopi, 2*n)

###############################################################################

# basis functions
def basis(l,m):
	return np.cos(theta)**m*np.sin(theta)**2*np.cos(2*l*theta)

# spectral Laplacian
def laplacian(l,m):
	return -4/(r*r)*(
		(l*l+m+1)*np.cos(2*l*theta) +
		l*(3.0-2*(m+2)*np.sin(theta)**2)*np.sin(2*l*theta)/np.sin(2*theta)
		)*np.cos(theta)**m

###############################################################################

# construct Laplacian operator matrix
B = np.zeros([n,n])
L = np.zeros([n,n])

# radial equation separates from phi!
for k in range(n):
	B[:,k] = basis(k,m)
	L[:,k] = laplacian(k,m)

# solve for radial eigenmodes
w,v = eig(L,B)

# index by eigenvalue
idx = np.argsort(-w.real)

'''
# output eigenvalues and residuals
for i in idx:
	print(-w[i].real,w[i].imag, norm(np.matmul(L,v[:,i]) - w[i]*np.matmul(B,v[:,i])))
'''

# 2D mesh LUTs
gamma = -w[idx[l]].real
f = np.matmul(B,v[:,idx[l]].real)
Phi,F = np.meshgrid(phi, f)

# symmetric color scale
z = np.max(np.abs(f))
levels = np.linspace(-z,z,100)

###############################################################################
# make a figure
###############################################################################

import matplotlib.pyplot as plt

# configure figure size
fig = plt.figure(figsize=(6, 4), frameon=False)

# set projection to polar coordinates
plt.subplot(111, polar=True)

# plot f(x) in plasma colormap, opaque
plt.contourf(phi, rho, F*np.cos(m*Phi), cmap='jet', extent=[0,twopi,0,r], levels=levels)
#plt.colorbar()

# manually set plot annotations
plt.title(r"$\gamma = %.6f$" % np.sqrt(gamma))
plt.xticks([])
plt.yticks([])

# tighten the whitespace (optional)
plt.tight_layout()

plt.show()
