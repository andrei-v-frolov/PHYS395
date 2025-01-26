#!/usr/bin/env python
# effects of collocation grid choice on expansion accuracy

#######################################################################

import numpy as np
from numpy.linalg import solve, svd, norm
from numpy.polynomial.legendre import legvander, legval, legder, legroots
from numpy.polynomial.chebyshev import chebvander, chebval, chebder, chebroots
from numpy.polynomial.polynomial import polyvander, polyval, polyder, polyroots

#######################################################################

# number of coefficients
n = 31

# allocate space for various collocation grids
grid = np.empty([6,n])

# uniform grid, endpoints INCLUDED vs EXCLUDED
dx = 2.0/n

grid[0] = np.linspace(-1.0, 1.0, n)
grid[3] = np.linspace(-1.0+dx/2, 1.0-dx/2, n)

# Legender grid, endpoints INCLUDED vs EXCLUDED
c = np.zeros(n+1); c[-1] = 1.0

grid[1] = np.concatenate([[-1.0], legroots(legder(c[1:])), [1.0]])
grid[4] = legroots(c)

# Chebyshev grid, endpoints INCLUDED vs EXCLUDED
dt = np.pi/n

grid[2] = np.cos(np.linspace(np.pi, 0.0, n))
grid[5] = np.cos(np.linspace(np.pi-dt/2, dt/2, n))

#######################################################################

# Legendre, Chebyshev, and monomial bases
P = legvander(grid, n-1)
T = chebvander(grid, n-1)
Q = polyvander(grid, n-1)

# compute condition number of a matrix
def condition(M):
	S = svd(M, compute_uv=False)
	return S[0]/S[-1]

# summary of collocation grids conditioning
print("Condition numbers for monomial basis:",  [condition(Q[i]) for i in range(6)])
print("Condition numbers for Legendre basis:",  [condition(P[i]) for i in range(6)])
print("Condition numbers for Chebyshev basis:", [condition(T[i]) for i in range(6)])

#######################################################################

# test function
f = np.exp(-grid**2*4.5)
#f = 1.0/(1.0 + 10.0*grid**2)

# expansion coefficients
c = solve(P,f)

# number of test points
pts = 1<<13

# evaluation grid
x = np.linspace(-1.0, 1.0, pts)
f = np.exp(-x*x*4.5)
#f = 1.0/(1.0 + 10.0*x*x)

# evaluate polynomial expansion and residual
g = legval(x,c.T)

# residual error metrics
print("Maximal residual:", np.max(np.abs(g-f),axis=1))
print("RMS error:", norm(g-f,axis=1)/pts)

#######################################################################

import matplotlib.pyplot as plt

fig = plt.figure(figsize=(8,3))

# label condition number on secondary y axis
ax = fig.gca(); con = ax.secondary_yaxis('right')

# closed and open collocation grids
for i,c in enumerate(['blue', 'orange', 'red']*2):
	plt.plot(grid[i], np.repeat(i,n), 'o', color = 'tab:'+c)
	#plt.plot(np.arccos(grid[i])/np.pi, np.repeat(i,n), 'o', color = 'tab:'+c)

# separator between closed and open grids
plt.axhline(2.5, color='black', linestyle='-', linewidth=0.5)

# text and condition number labels
plt.yticks(range(6), ['uniform', 'Legendre', 'Chebyshev']*2)
con.set_yticks(range(6), [f'{condition(T[i]):.3g}' for i in range(6)])

# restrict x axis range
plt.xlim([-1.0,1.0])
#plt.xlim([0.0,1.0])

plt.show()
