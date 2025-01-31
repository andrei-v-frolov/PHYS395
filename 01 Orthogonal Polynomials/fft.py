#!/usr/bin/env python
# compute (complex) fast Fourier transform of a function

#######################################################################

import numpy as np
from numpy.fft import fft, ifft, fftfreq

# number of grid points
n = 256

# uniform evaluation grid - EXCLUDING the right endpoint
l = 1.0; dx = 2.0*l/n
x = np.linspace(-l, l-dx, n)

# wavenumber grid using helper function
nn = n//2; dk = np.pi/l
k = (2.0*np.pi) * fftfreq(n, dx)

# or, explicitly compute wavenumber grid
#k = np.fromiter(((i-n if i>nn else i)*dk for i in range(0,n)), 'double')

# note that helper returns NEGATIVE of Nyquist frequency 
print(f'Nyquist frequency is {k[n//2]}')

#######################################################################

# GENERIC test function (and its derivative)
f = np.exp(-x*x*4.5); df = -9.0*x*f

# PERIODIC test function (and its derivative)
#y = np.sin(np.pi*x/(2*l)); z = np.cos(np.pi*x/(2*l)) 
#f = np.exp(-y*y*4.5); df = -(9.0*np.pi/(2*l))*y*z*f

# NYQUIST mode function (and its derivative)
#f = np.cos(np.pi*nn*x); ddf = -(np.pi*nn)**2*f

#######################################################################

# compute (complex) Fourier transform
F = fft(f); g = ifft(F)

# compute derivatives
dg = ifft(1.0j*k*F)
ddg = ifft(-k*k*F)

#######################################################################

import matplotlib.pyplot as plt

#plt.plot(x, f)
plt.plot(x, np.real(dg-df), "r-")
plt.plot(x, np.imag(dg-df), "b-")

# restrict x axis range
plt.xlim([-l,l])

# show in interactive console
plt.show()