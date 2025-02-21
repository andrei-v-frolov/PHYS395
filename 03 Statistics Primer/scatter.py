#!/usr/bin/env python
# scatter plot of bivariate PDF

#######################################################################

# standard stream access
from sys import stdin

# numerical libraries
import numpy as np

######################################################################## read in data from stdin

# read in data from stdin
data = np.loadtxt(stdin)
n,columns = data.shape

# sanity check on data format
assert columns > 1, "Expecting at least two columns of data..."

#######################################################################

import matplotlib.pyplot as plt

# extract data to plot
x = data[:,0]   # Apple
y = data[:,11]  # Starbucks
z = data[:,8]   # Netflix

# scatter plot, with (optional) color coding
plt.scatter(x, y, c=np.arange(n), cmap='plasma', marker=".", s=3)

# binned histograms, with log(pdf) used for color
#plt.hist2d(x, y, bins=30, norm='log', cmap='Blues')
#plt.hexbin(x, y, gridsize=30, bins='log', cmap='Blues')

# set square aspect if desired
#plt.gca().set_aspect('equal')

plt.colorbar()
plt.show()
