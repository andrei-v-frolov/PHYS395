#!/usr/bin/env python
# compute fast Chebyshev transform of a function

import numpy as np
from scipy.fft import dct, idct, idst

# number of points
n = 1024

# evaluation grids
dt = np.pi/n; theta = np.linspace(np.pi-dt/2, dt/2, n)
l = 1.0; x = l * np.cos(theta); w = l * np.sin(theta)
k = np.linspace(0.0, n, n, endpoint=False)

# function to approximate (and its derivative)
f = np.exp(-x*x*4.5); df = -9.0*x*f

# compute Type-II discrete cosine transform
F = dct(f); g = idct(F)

# compute Chebyshev derivative using Type-II DST
dg = -idst(np.roll(k*F,-1))/w

#######################################################################

import matplotlib.pyplot as plt

# plot function residual
#plt.plot(x,f,'b.')
#plt.plot(x,g,'r-')

# plot derivative residual
#plt.plot(x,df,'b.')
plt.plot(x,df-dg,'r-')

plt.show()
