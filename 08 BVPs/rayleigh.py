#!/usr/bin/env python
# find energy eigenstates directly and using Rayleigh iteration

#######################################################################

import numpy as np
from numpy.linalg import solve

#######################################################################

# number of modes and compactification scale
n = 201; l = 1.0

# Chebyshev collocation grid
dt = np.pi/n; t = np.linspace(np.pi-dt/2.0, dt/2.0, n)

# construct Laplacian operator matrix
B = np.zeros([n,n])
D = np.zeros([n,n])

for k in range(0,n):
	B[k] = np.cos(k*t)
	D[k] = -k * (k*np.sin(t)*np.cos(k*t) + 2*np.cos(t)*np.sin(k*t)) * np.sin(t)**3/l**2

L = solve(B,D).T

#######################################################################

from numpy.linalg import eig

# rational collocation grid
x = l/np.tan(t)

# potential of a quantum oscillator
def V(x):
	return x*x/2.0

# Hamiltonian of a quantum oscillator
H = -L/2.0 + np.diag(V(x))

# eigenvalues and eigenfunctions
sigma,psi = eig(H); idx = np.argsort(sigma)

#######################################################################

# scalar product weights
w = l/np.sin(t)**2

# Rayleigh quotient
def rayleigh(u):
	return np.sum(w*u*(H@u))/np.sum(w*u*u)

# normalize wavefunction to unit norm and definite sign
def normalize(u):
	return u/np.sqrt(np.sum(w*u*u))*np.sign(np.sum(w*(1.0+x)*u))

# intial guess
phi = x*np.exp(-x*x/1.0)

# iteration history
history = np.zeros([16,n])

# Rayleigh iteration
for i in range(0,16):
	phi = normalize(phi); history[i] = phi
	mu = rayleigh(phi); print(mu)
	phi = solve(H-mu*np.eye(n), phi)

#######################################################################

import matplotlib.pyplot as plt
import matplotlib.colors as clr

'''
# Hamiltonian operator matrix
plt.imshow(H, origin='lower', cmap='cividis', norm=clr.SymLogNorm(100.0, vmin=-1.0e3,vmax=1.0e3))
plt.colorbar()
'''

'''
# plot energy eigenstates
fig, ax = plt.subplots(1, 2, sharey=True)
for k in range(n):
	ax[0].plot(x, sigma[idx[k]]+2.0*psi[:,idx[k]], color='black')
ax[0].plot(x, V(x), color='red', linewidth=3)
ax[0].set_xlim([-5,5])

ax[1].plot(range(n), sigma[idx])
ax[1].set_xlim([0,(n-1)//16])

plt.ylim([0,sigma[idx[(n-1)//16]]+0.5])
plt.tight_layout()
'''

# Rayleigh iterations
plt.plot(x,history.T)
plt.xlim([-5,5])

plt.show()
