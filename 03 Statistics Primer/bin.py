#!/usr/bin/env python
# bin IID random samples to estimate PDF

#######################################################################

# operating system functions (for standard stream access)
from sys import stdin

# numerical libraries
import numpy as np

#######################################################################

# read in data from stdin
data = np.loadtxt(stdin)

# sanity check on data format
assert len(data.shape) == 1, "Expecting single column data..."

# data bounds
a = np.min(data)
b = np.max(data)

# number of data bins
n = 32; dx = (b-a)/n; x = np.linspace(a+dx/2.0, b-dx/2.0, n)

# compute PDF histogram, normalizing density
hist,edge = np.histogram(data, bins=n, range=(a,b), density=True)

#######################################################################

import matplotlib.pyplot as plt

plt.bar(x, hist, width=0.8*dx, align='center')

#######################################################################

'''
from scipy.interpolate import PchipInterpolator as interpolate

# smooth interpolation with monotonic spline
pdf = interpolate(x, hist)
grid = np.linspace(a, b, 1024)
plt.fill_between(grid, pdf(grid))
'''

plt.ylim(bottom=0.0)

plt.show()
