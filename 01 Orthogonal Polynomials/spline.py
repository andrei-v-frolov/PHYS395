#!/usr/bin/env python
# compute cubic spline approximation

#######################################################################

import numpy as np

# interval to sample on
n = 8; a = -np.pi; b = np.pi

x = np.linspace(a, b, n)
y = np.sin(x)

#######################################################################

from scipy.interpolate import interp1d, CubicSpline

# legacy interpolation interface
spline = interp1d(x, y, kind='cubic')

# new style drivers are more flexible
#spline = CubicSpline(x, y, bc_type='periodic')

x = np.linspace(a, b, 1024)
y = np.sin(x); f = spline(x)

#######################################################################

#from pltconfig import *
import matplotlib.pyplot as plt

plt.figure()

plt.plot(x, f, '-')
plt.axhline(0.0, color='black', linestyle=':', zorder=0)
plt.fill_between(x, y, 0.0, color='green', alpha=0.2)

plt.plot(x, f-y, 'r-')

plt.xlim([a,b])

# show in interactive console
plt.show()

# render the figure to a file
#plt.savefig("output.pdf", bbox_inches='tight', pad_inches=0.02, transparent=True)
