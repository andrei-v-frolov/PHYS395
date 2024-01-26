#!/usr/bin/env python
# compute fast Fourier transform of a function

import numpy as np
from scipy.fft import fft, ifft

# number of points
n = 32

# evaluation grid - uniform
l = 1.0; dx = 2.0*l/n
x = np.linspace(-1.0, 1.0-dx, n)

# function to approximate
f = np.exp(-x*x*4.5)

# compute (complex) Fourier transform
F = fft(f); g = ifft(F)

#######################################################################

import matplotlib.pyplot as plt

# plot residual
plt.plot(x,f,'b.')
plt.plot(x,np.real(g),'r-')

plt.show()
