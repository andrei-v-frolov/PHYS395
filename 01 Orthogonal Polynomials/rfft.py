#!/usr/bin/env python
# compute (real) fast Fourier transform of a function

import numpy as np
from scipy.fft import rfft, irfft, rfftfreq

# number of points
n = 1024

# evaluation grid - uniform
l = 1.0; dx = 2.0*l/n
x = np.linspace(-l, l-dx, n)
k = (2.0*np.pi) * rfftfreq(n, dx)

# function to approximate (and its derivative)
f = np.exp(-x*x*4.5); df = -9.0*x*f

# compute (real) Fourier transform
F = rfft(f); g = irfft(F); dg = irfft(1.0j*k*F)

#######################################################################

import matplotlib.pyplot as plt

# plot function residual
#plt.plot(x,f,'b.')
#plt.plot(x,g,'r-')

# plot derivative residual
#plt.plot(x,df,'b.')
plt.plot(x,df-dg,'r-')

plt.show()
