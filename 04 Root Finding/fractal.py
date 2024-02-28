#!/usr/bin/env python
# basin of convergence for Newton's method (fractal for z^3=1)

#######################################################################

import numpy as np

#######################################################################

# resolution and interval to sample
n = 1024; l = 3.0

# number of iterations needed to converge to given tolerance
def capture(z, iterations = 50, epsilon = 1.0e-6):
	for i in range(0,iterations):
		f = z**3-1.0; residual = abs(f)
		if (residual < epsilon): break
		z -= f/(3.0*z*z)
	# since convergence is quadratic, we can extrapolate to floats
	return i - np.log2(np.log(min(residual,epsilon))/np.log(epsilon))

# grid of starting values
x = np.linspace(-l, l, n)
X,Y = np.meshgrid(x,x)

# evaluate convergence for all starting values
q = np.vectorize(capture)(X + 1.0j*Y)

#######################################################################

'''
# fast decimator using CIC filter (for oversampled rendering)
# https://en.wikipedia.org/wiki/Cascaded_integrator–comb_filter
q = np.diff(np.cumsum(q, axis=0)[::4,:], axis=0)/4.0
q = np.diff(np.cumsum(q, axis=1)[:,::4], axis=1)/4.0
'''

#######################################################################

import matplotlib.pyplot as plt

plt.imshow(q, origin='lower', extent=[-l,l,-l,l], interpolation='none', norm='log', cmap='twilight')
plt.colorbar()

plt.show()