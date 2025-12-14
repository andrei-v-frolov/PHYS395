#!/usr/bin/env python
# compute (real) fast Fourier transform of a function of two variables

#######################################################################

import numpy as np
from numpy.fft import rfft2, irfft2, rfftfreq, fftfreq

#######################################################################

# number of grid points
nx = 256; ny = 512

# uniform evaluation grid - EXCLUDING the right endpoint
lx = 1.0; dx = 2.0*lx/nx; x = np.linspace(-lx, lx-dx, nx)
ly = 2.0; dy = 2.0*ly/ny; y = np.linspace(-ly, ly-dy, ny)

# wavenumber grids using helper functions
kx = (2.0*np.pi) * rfftfreq(nx, dx)
ky = (2.0*np.pi) * fftfreq(ny, dy)

# 2D grid iterators
X,Y = np.meshgrid(x,y)
Kx,Ky = np.meshgrid(kx,ky)

# Laplacian operator
K2 = Kx*Kx + Ky*Ky

#######################################################################

# generic test function (and its derivatives)
f = np.exp(-X*X*45-Y*Y*15)
dxf = -90.0*X*f; dyf = -30.0*Y*f
ddf = 60.0*(135*X*X + 15*Y*Y - 2.0)*f

#######################################################################

# compute (real) Fourier transform
F = rfft2(f); g = irfft2(F)

# compute derivatives
dxg = irfft2(1.0j*Kx*F)
dyg = irfft2(1.0j*Ky*F)
ddg = irfft2(-K2*F)

#######################################################################

import matplotlib.pyplot as plt

fig = plt.figure(); ax = fig.gca()

plt.imshow(ddg-ddf, extent=[-lx,lx,-ly,ly], origin='lower', cmap='bwr', interpolation='none')
plt.colorbar()

# show in interactive console
plt.show()
