#!/usr/bin/env python
# compute cubic spline approximation

#######################################################################

import numpy as np
from scipy.interpolate import CubicSpline

#######################################################################

# number of points
n = 20

# evaluation grid - uniform
dx = 2.0*np.pi/n
x = np.linspace(dx/2.0, 2.0*np.pi-dx/2.0, n)

# function to approximate
f = np.sin(x)

# cubic spline interpolation
spline = CubicSpline(x,f)

# exploiting extra information on BCs
#x = np.linspace(0.0, 2.0*np.pi, n); f = np.sin(x)
#spline = CubicSpline(x,f, bc_type='periodic')

#######################################################################

import matplotlib.pyplot as plt

# test grid
pts = 1000
x = np.linspace(0.0, 2.0*np.pi, pts)
f = np.sin(x)
y = spline(x)

# plot residual
#plt.plot(x,f,'b-')
plt.plot(x,y-f,'r-')

plt.show()
