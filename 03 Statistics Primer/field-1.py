#!/usr/bin/env python
# generate a Gaussian random field in 1D using convolution

#######################################################################

import numpy as np
from numpy.random import normal
from scipy.fft import rfft, irfft

#######################################################################

# number of points
n = 1024

# uniform evaluation grid
l = 1.0; dx = 2.0*l/n
x = np.linspace(-l, l-dx, n)

# white random field and window function
f = normal(size=n); w = np.exp(-x*x*450.0)

# window function can be normalized if you want
w /= np.sum(w)

# convolution can be done directly...
g = np.convolve(f,w, mode='same')

# ... or using (real) Fourier transform
F = rfft(f); W = rfft(w); q = irfft(W*F)

#######################################################################

import matplotlib.pyplot as plt

# plot random field
plt.plot(x,f,'b-')

# note different window center locations!
plt.plot(x,np.roll(g,-1),'k-')
plt.plot(x,np.roll(q,n//2),'r-')

plt.show()
