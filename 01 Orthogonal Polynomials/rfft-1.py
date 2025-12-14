#!/usr/bin/env python
# compute (real) fast Fourier transform of a function

#######################################################################

import numpy as np
from numpy.fft import rfft, irfft, rfftfreq

# number of grid points
n = 256

# uniform evaluation grid - EXCLUDING the right endpoint
l = 1.0; dx = 2.0*l/n
x = np.linspace(-l, l-dx, n)

# wavenumber grid using helper function
nn = n//2; dk = np.pi/l
k = (2.0*np.pi) * rfftfreq(n, dx)

# or, explicitly compute wavenumber grid
#k = np.arange(nn+1)*dk

# note that helper returns POSITIVE of Nyquist frequency 
print(f'Nyquist frequency is {k[n//2]}')

#######################################################################

# GENERIC test function (and its derivative)
f = np.exp(-x*x*4.5); df = -9.0*x*f

# PERIODIC test function (and its derivative)
#w = np.pi/(2*l); y = np.sin(w*x)/w
#f = np.exp(-y*y*4.5); df = -9.0*y*f*np.cos(w*x)

# NYQUIST mode function (and its derivative)
#f = np.cos(np.pi*nn*x); ddf = -(np.pi*nn)**2*f

#######################################################################

# compute (real) Fourier transform
F = rfft(f); g = irfft(F)

# compute derivatives
dg = irfft(1.0j*k*F)
ddg = irfft(-k*k*F)

#######################################################################

import matplotlib.pyplot as plt

#plt.plot(x, f)
plt.plot(x, dg-df, "r-")

# restrict x axis range
plt.xlim([-l,l])

# show in interactive console
plt.show()