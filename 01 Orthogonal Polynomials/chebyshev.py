#!/usr/bin/env python
# compute fast Chebyshev transform of a function

#######################################################################

import numpy as np
from scipy.fft import dct, idct, idst

# number of grid points
n = 256

# cosine evaluation grid - EXCLUDING endpoints
dt = np.pi/n; theta = np.linspace(np.pi-dt/2, dt/2, n)
l = 1.0; x = l*np.cos(theta); k = np.arange(n)/l

#######################################################################

# GENERIC test function (and its derivative)
f = np.exp(-x*x*4.5); df = -9.0*x*f

# PERIODIC test function (and its derivative)
#w = np.pi/(2*l); y = np.sin(w*x)/w
#f = np.exp(-y*y*4.5); df = -9.0*y*f*np.cos(w*x)

#######################################################################

# compute Type-II discrete cosine transform
F = dct(f); g = idct(F)

# compute Chebyshev derivative using Type-II DST
dg = -idst(np.roll(k*F,-1))/np.sin(theta)

#######################################################################

import matplotlib.pyplot as plt

#plt.plot(x, f)
plt.plot(x, dg-df, "r-")

# restrict x axis range
plt.xlim([-l,l])

# show in interactive console
plt.show()