#!/usr/bin/env python
# bin IID random samples to estimate PDF

#######################################################################

# standard input access
from sys import stdin

# numerical libraries
import numpy as np

#######################################################################

# read in data from stdin
data = np.loadtxt(stdin)

# sanity check on data format
assert len(data.shape) == 1, "Expecting single column data..."

# number of bins (could be estimated from data, e.g. n='fd')
n = 32

# data bounds
a = data.min()
b = data.max()

# compute PDF histogram, normalizing density
pdf,edges = np.histogram(data, bins=n, range=(a,b), density=True)

# bin centers and widths
x = (edges[1:]+edges[:-1])/2.0; dx = np.diff(edges)

#######################################################################

import matplotlib.pyplot as plt

#plt.stairs(pdf, edges)
plt.bar(x, pdf, width=0.8*dx, align='center')
#plt.stairs(np.cumsum(pdf*dx), edges)

'''
# smooth interpolation with monotonic spline
from scipy.interpolate import PchipInterpolator as interp
grid = np.linspace(a, b, 1024)
plt.fill_between(grid, interp(x, pdf)(grid))
'''

# restrict axis range
plt.xlim([a,b]); plt.ylim(bottom=0.0)

# show in interactive console
plt.show()
