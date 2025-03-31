#!/usr/bin/env python
# time finite difference Laplacian stencils precompiled using Numba

#######################################################################

import numpy as np
from numba import jit
from time import perf_counter as now

#######################################################################

n = 1000000

# test data for stencil application
x = np.zeros(n); x[3] = 1.0

# 1D laplacian stencil via traditional loop
@jit(nogil=True,cache=True)
def laplacian_1d(x):
	n, = x.shape; y = np.zeros_like(x)
	for i in range(1,n-1):
		y[i] = x[i-1] - 2.0*x[i] + x[i+1]
	return y

t1 = now(); y = laplacian_1d(x); t2 = now()
#print("numba compile (1D):    ", y[0:7], t2-t1)

t1 = now(); y = laplacian_1d(x); t2 = now()
print("numba execute (1D):    ", y[0:7], t2-t1)

#######################################################################

n = 1000

# test data for stencil application
x = np.zeros([n,n]); x[3,3] = 1.0

# 2D laplacian stencil via traditional loop
@jit(nogil=True,cache=True)
def laplacian_2d(x):
	m,n = x.shape; y = np.zeros_like(x)
	for i in range(1,m-1):
		for j in range(1,n-1):
			y[i,j] = x[i-1,j] + x[i,j-1] - 4.0*x[i,j] + x[i,j+1] + x[i+1,j]
	return y

t1 = now(); y = laplacian_2d(x); t2 = now()
#print("numba compile (2D):    ", y[0:7,0:7], t2-t1, sep='\n')

t1 = now(); y = laplacian_2d(x); t2 = now()
print("numba execute (2D):    ", y[0:7,0:7], t2-t1, sep='\n')
