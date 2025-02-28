#!/usr/bin/env python
# demo of built-in SciPy optimization methods

#######################################################################

from scipy.optimize import minimize

#######################################################################

# tilt parameter
mu = 0.1

# call history
history = []

# tilted Mexican hat (testing convergence in winding valley)
def f(v):
	global history; history.append(v)
	x,y = v; return (x*x+y*y - 1.0)**2/4.0 - mu*x

def df(v):
	x,y = v; return [(x*x+y*y - 1.0)*x - mu, (x*x+y*y - 1.0)*y]

#######################################################################

# starting point
start = [0.0,2.0]

# if not supplied, gradient is evaluated by finite difference
result = minimize(f, start); print("SciPy default", result)
result = minimize(f, start, method='BFGS', jac=df); print("BFGS", result)
result = minimize(f, start, method='L-BFGS-B', jac=df); print("L-BFGS-B", result)
result = minimize(f, start, method='CG', jac=df); print("Conjugate gradient", result)
result = minimize(f, start, method='Powell'); print("Powell's method", result)
result = minimize(f, start, method='Nelder-Mead'); print("Nelder-Mead method", result)

#######################################################################

import numpy as np
import matplotlib.pyplot as plt

# domain to be plotted
l = 1.5; pts = 1024; dl = l/pts

# function being minimized
grid = np.linspace(-l+dl,l-dl,pts)
X,Y = np.meshgrid(grid,grid)
F = np.vectorize(lambda x,y: f([x,y]))(X,Y)

# scale function for plotting
F = F/(mu*mu+F*F)**0.375

# re-run specific optimizer to plot
history = []
minimize(f, start, method='BFGS', jac=df)
history = np.array(history)

# plot optimization history
plt.plot(history[:,0], history[:,1], "r.-")
#plt.imshow(F, extent=[-l,l,-l,l], cmap='Blues')
plt.contourf(F, extent=[-l,l,-l,l], levels=13, cmap='Blues')
plt.colorbar()

# make sure aspect is 1:1
plt.gca().set_aspect('equal')
plt.xlim([-l,l]); plt.ylim([-l,l])

plt.show()
