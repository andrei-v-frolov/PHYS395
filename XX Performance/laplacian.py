#!/usr/bin/env python
# time evaluating finite difference Laplacian in different ways

#######################################################################

import numpy as np
import scipy.signal as sig
import scipy.ndimage as img

from time import perf_counter as now

#######################################################################

n = 1000000

# test data for stencil application
x = np.zeros(n); x[3] = 1.0

# second order Laplacian stencil in 1D
stencil = np.array([1,-2,1])

t1 = now(); y = np.convolve(x,stencil, mode='same'); t2 = now()
print("numpy.convolve:        ", y[0:7], t2-t1)

t1 = now(); y = x[:-2]-2*x[1:-1]+x[2:]; t2 = now()
print("numpy slice arrays:    ", y[0:7], t2-t1)

t1 = now(); y = np.roll(x,-1)-2*x+np.roll(x,1); t2 = now()
print("numpy roll arrays:     ", y[0:7], t2-t1)

t1 = now(); y = sig.convolve(x,stencil, mode='same', method='direct'); t2 = now()
print("scipy.signal.convolve: ", y[0:7], t2-t1)

t1 = now(); y = img.convolve(x,stencil, mode='wrap'); t2 = now()
print("scipy.ndimage.convolve:", y[0:7], t2-t1)

#######################################################################

n = 1000

# test data for stencil application
x = np.zeros([n,n]); x[3,3] = 1.0

# second order Laplacian stencil in 2D
stencil = np.array([[0,1,0],[1,-4,1],[0,1,0]])

# second order isotropic Laplacian stencil in 2D
#stencil = np.array([[1,4,1],[4,-20,4],[1,4,1]])/6.0

t1 = now(); y = x[1:-1,:-2]+x[:-2,1:-1]-4*x[1:-1,1:-1]+x[2:,1:-1]+x[1:-1,2:]; t2 = now()
print("numpy slice arrays:", y[0:7,0:7], t2-t1, sep='\n')

t1 = now(); y = img.convolve(x,stencil, mode='wrap'); t2 = now()
print("scipy.ndimage.convolve:", y[0:7,0:7], t2-t1, sep='\n')

t1 = now(); y = sig.convolve(x,stencil, mode='same', method='direct'); t2 = now()
print("scipy.signal.convolve:", y[0:7,0:7], t2-t1, sep='\n')

#######################################################################
