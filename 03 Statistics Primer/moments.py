#!/usr/bin/env python
# estimate moments and L-moments of sampled data

#######################################################################

# operating system functions (for standard stream access)
from sys import stdin

# numerical libraries
import numpy as np
from numpy.polynomial.legendre import legvander

#######################################################################

# read in data from stdin
data = np.loadtxt(stdin)

# sanity check on data format
assert len(data.shape) == 1, "Expecting single column data..."

# average-subtracted data for estimating central moments
n = len(data); avg = np.sum(data)/n; x = data-avg

# sorted data and Legendre polynomials for estimating L-moments
y = np.sort(data); P = legvander(np.linspace(-1.0, 1.0, n), 5)

# output computed moments, same order per line
for k in range(1,5):
	m = np.sum(x**k)/n
	l = np.sum(y*P[:,k-1])/n
	print(m, l)
