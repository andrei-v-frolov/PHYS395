#!/usr/bin/env python
# convergence of Newton's method (fractal for z^3=1)

#######################################################################

# numerical libraries
import numpy as np

# just-in-time compiler (comment out if unavailable)
from numba import vectorize, float64, complex128

#######################################################################

# resolution, oversampling, and interval to sample
n = 1024; os = 4; l = 3.0
n = n*os + (2 if os > 1 else 0)

# grid of starting values
re = np.linspace(-l, l, n)
im = np.linspace(-l, l, n)
a,b = np.meshgrid(re,im)

# form a complex-valued array
z = a + b*1j

#######################################################################

'''
# simplest way to see the structure
for i in range(8):
	z -= (z**3-1)/(3.0*z**2)

# not the prettiest, but servicable
q = np.abs(z**3-1) + 1.0e-6
'''

#######################################################################

# max iterations and tolerance
iterations = 50; epsilon = 1.0e-6

# number of iterations needed to converge
@vectorize([float64(complex128)])
def capture(z):
	for i in range(iterations):
		f = z**3-1.0; residual = abs(f)
		if (residual < epsilon): break
		z -= f/(3.0*z*z)
	# since convergence is quadratic, we can extrapolate it to floats
	return i - np.log2(np.log(min(residual,epsilon))/np.log(epsilon))

# evaluate convergence for all starting values
#q = np.vectorize(capture)(z)
q = capture(z)

#######################################################################

# fast decimator using CIC filter (for oversampled rendering)
# https://en.wikipedia.org/wiki/Cascaded_integrator–comb_filter
q = np.diff(np.cumsum(q, axis=0)[::os,:], axis=0)/os
q = np.diff(np.cumsum(q, axis=1)[:,::os], axis=1)/os

#######################################################################

import matplotlib.pyplot as plt

plt.imshow(q, origin='lower', extent=[-l,l,-l,l], norm='log', cmap='twilight', interpolation='none')

# set square aspect if desired
#plt.gca().set_aspect('equal')
#plt.xlim([-l,l]); plt.ylim([-l,l])

plt.colorbar()
plt.show()

