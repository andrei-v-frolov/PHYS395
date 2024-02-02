#!/usr/bin/env python
# draw a number of IID samples from a distribution
# run as: python draw.py PDF [dof:k] [samples:n] [seed:x]

#######################################################################

# default argument values
pdf = "default"; initial = None; n = 32; dof = 7

#######################################################################

import numpy as np
from numpy.random import seed, rand, uniform, normal

# seed random sequence
seed(initial)

# draw samples from specified distribution
match pdf:
	case "uniform":
		x = uniform(size=n)
	case "normal":
		x = normal(size=n)
	case "sum":
		x = np.sum(uniform(size=[n,dof]), 1)
	case "max":
		x = np.max(normal(size=[n,dof]), 1)
	case "chi2":
		x = np.sum(normal(size=[n,dof])**2, 1)
	case _:
		x = rand(n)

# output samples, one per line, in scientific notation
for value in x:
	print("%.16e" % value)
