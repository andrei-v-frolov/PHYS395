#!/usr/bin/env python
# scatter plot of bivariate PDF, with a few variations

#######################################################################

# operating system functions (for standard stream access)
from sys import stdin

# numerical libraries
import numpy as np

#######################################################################

# read in data from stdin
data = np.loadtxt(stdin)

# sanity check on data format
assert len(data.shape) > 1, "Expecting at least two columns of data..."

#######################################################################

import matplotlib.pyplot as plt

x = data[:,0]	# Apple
y = data[:,11]	# Starbucks
z = data[:,8]	# Netflix

t = np.linspace(0.0, 1.0, len(data))

# scatter plot, with (optional) color coding
#plt.scatter(x, y, c=t, cmap='plasma', marker=".")
plt.scatter(x, y, c=z, cmap='plasma', marker=".")

# binned histograms, with log(pdf) used for color
#plt.hist2d(x, y, bins=30, norm='log', cmap='Blues')
#plt.hexbin(x, y, gridsize=30, bins='log', cmap='Blues')

plt.colorbar()
plt.show()
