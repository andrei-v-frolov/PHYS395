#!/usr/bin/env python
# compute (complex) fast Fourier transform of a function

import numpy as np
from scipy.fft import fft, ifft

# number of points
n = 1024

# evaluation grid - uniform
l = 1.0; dx = 2.0*l/n
x = np.linspace(-l, l-dx, n)

# function to approximate (and its derivative)
f = np.exp(-x*x*4.5); df = -9.0*x*f

# compute (complex) Fourier transform
F = fft(f); g = ifft(F)

# compute frequency grid
nn = n//2; dk = np.pi/l

k = np.zeros(n)
for i in range(0,nn+1):
	k[i] = i*dk
for i in range(nn+1,n):
	k[i] = (i-n)*dk

# compute derivative
dg = ifft(1.0j*k*F)

#######################################################################

import matplotlib.pyplot as plt

# plot function residual
#plt.plot(x,f,'b.')
#plt.plot(x,np.real(g),'r-')

# plot derivative residual
#plt.plot(x,df,'b.')
plt.plot(x,df-np.real(dg),'r-')

plt.show()
