#!/usr/bin/env python
# draw a specified number of IID samples from a distribution
# run as: python draw.py PDF [dof:k] [samples:n] [seed:x]

#######################################################################

from sys import argv

# default argument values
pdf = "default"; initial = None; n = 32; dof = 7

# select distribution (if specified)
if len(argv) > 1:
	pdf = argv[1]

# parse parameters (if supplied)
for p in argv[2:]:
	match tuple(p.split(":")):
		case ("dof",s):
			dof = int(s)
		case ("samples",s):
			n = int(s)
		case ("seed",s):
			initial = int(s)
		case _:
			...

#######################################################################

import numpy as np
from erfinv import erfinv
from numpy.random import seed, rand, uniform, normal

# seed random sequence
seed(initial)

# draw samples from specified distribution
match pdf:
	# built-in generators
	case "uniform":
		x = uniform(size=n)
	case "normal":
		x = normal(size=n)
	# aggregate distributions
	case "sum":
		x = np.sum(uniform(size=[n,dof]), 1)
	case "max":
		x = np.max(normal(size=[n,dof]), 1)
	case "chi2":
		x = np.sum(normal(size=[n,dof])**2, 1)
	# shaping uniform PDF
	case "erfinv":
		x = np.sqrt(2.0) * np.vectorize(erfinv)(2.0*rand(n)-1.0)
	case "box-muller":
		u = rand(n); v = rand(n)
		x = np.sqrt(-2.0*np.log(u)) * np.cos(2.0*np.pi*v)
	# default random generator
	case _:
		x = rand(n)

# output samples, one per line, in scientific notation
for value in x:
	print("%.16e" % value)
