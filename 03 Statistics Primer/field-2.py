#!/usr/bin/env python
# generate a Gaussian random field in 2D using convolution

#######################################################################

import numpy as np
from numpy.random import normal
from scipy.fft import rfft2, irfft2

#######################################################################

# number of points in [n,n] grid
n = 1024

# uniform evaluation grid
l = 1.0; dx = 2.0*l/n
x = np.linspace(-l, l-dx, n)

# 2D grid iterators
X,Y = np.meshgrid(x,x)


# white random field and window function
f = normal(size=[n,n]); w = np.exp(-(X*X+Y*Y)*450.0)

# window function can be normalized if you want
w /= np.sum(w)

# convolution using (real) Fourier transform
F = rfft2(f); W = rfft2(w); g = irfft2(W*F)

#######################################################################

import matplotlib.pyplot as plt

plt.imshow(g, cmap='bwr', origin='lower')

plt.show()
