#!/usr/bin/env python
# convergence of Newton's method (fractal for z^3=1)

#######################################################################

# numerical libraries
import numpy as np

# just-in-time compiler (comment out if unavailable)
from numba import vectorize, float64, complex128

#######################################################################

# resolution, oversampling, and interval to sample
n = 1024; os = 1; l = 3.0; dx = l/(n*os)

# grid of starting values
re = np.linspace(-l+dx, l-dx, n*os)
im = np.linspace(-l+dx, l-dx, n*os)
a,b = np.meshgrid(re,im)

# form a complex-valued array
z = a + b*1j

#######################################################################
# Newton's fractal (convergence of Newton's method solving z^3=1)
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
# Mandelbrot fractal (iterating z -> z^2+c until divergence)
#######################################################################

'''
# max iterations and escape radius
iterations = 100; radius = 10.0

# number of iterations needed to escape
@vectorize([float64(complex128)])
def escape(c):
	z = 0j
	for i in range(iterations):
		z = z*z + c; r = abs(z)
		if (r > radius): break
	return i + np.sqrt(0.5) - np.log2(np.log(r))

# evaluate convergence for all starting values
#q = np.vectorize(escape)(z)
q = escape(z)

# replace NaNs with 0.0 for averaging
#q = np.nan_to_num(q)
'''

#######################################################################

'''
# fast decimator using CIC filter (for oversampled rendering)
# https://en.wikipedia.org/wiki/Cascaded_integrator–comb_filter
q = np.pad(q, (1,0))
q = np.diff(np.cumsum(q, axis=0)[::os,:], axis=0)/os
q = np.diff(np.cumsum(q, axis=1)[:,::os], axis=1)/os
'''

'''
# alternatively, median filter could be used instead
from scipy.ndimage import median_filter as median
q = median(q, os)[os//2::os,os//2::os]
'''

#######################################################################

import matplotlib.pyplot as plt

plt.imshow(q, origin='lower', extent=[-l,l,-l,l], norm='log', cmap='twilight', interpolation='none')

# set square aspect if desired
#plt.gca().set_aspect('equal')
#plt.xlim([-l,l]); plt.ylim([-l,l])

plt.colorbar()
plt.show()

