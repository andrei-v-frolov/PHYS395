#!/usr/bin/env python
# visualize vector field with streamlines

#######################################################################

import numpy as np

#######################################################################

# grid resolution and box size
n = 1024; l = 2.0

# uniform spatial grid
dx = 2.0*l/n; x = np.linspace(-l,l,n)

# electric field of a dipole
X,Y = np.meshgrid(x,x)
X1 = X-1.0; X2 = X+1.0
R1 = (X1*X1+Y*Y)**1.5
R2 = (X2*X2+Y*Y)**1.5
Ex = X1/R1 - X2/R2
Ey = Y/R1 - Y/R2

#######################################################################

import matplotlib.pyplot as plt

plt.gca().set_aspect('equal')
plt.streamplot(X,Y,Ex,Ey)

plt.show()
