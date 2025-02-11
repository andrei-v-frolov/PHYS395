#!/usr/bin/env python
# draw a specified number of IID samples from a distribution
# run as: python draw.py [distribution] [dof:k] [n:samples] [seed:value]

#######################################################################

# access command line arguments
from sys import argv

# default argument values
pdf = "default"; dof = 7
n = 1024; init = None

# parse command line arguments (if supplied)
for arg in argv[1:]:
	match tuple(arg.split(":")):
		case s if len(s) == 1:
			pdf, = s
		case ("dof",s):
			dof = int(s)
		case ("n",s):
			n = int(s)
		case ("seed",s):
			init = int(s)
		case _:
			...

#######################################################################

# bulk random number generators from NumPy
import numpy as np
from erfinv import erfinv
from numpy.random import seed, rand, uniform, normal

# seed random sequence
seed(init)

# draw samples from specified distribution
match pdf:
	# built-in generators
	case "uniform":
		x = uniform(size=n)
	case "normal":
		x = normal(size=n)
	# aggregate distributions
	case "sum":
		x = np.sum(uniform(size=[n,dof]), axis=1)
	case "max":
		x = np.max(normal(size=[n,dof]), axis=1)
	case "chi2":
		x = np.sum(normal(size=[n,dof])**2, axis=1)
	# shaping uniform PDF
	case "cauchy":
		x = np.tan(np.pi*(rand(n)-0.5))
	case "erfinv":
		x = np.sqrt(2.0) * erfinv(2.0*rand(n)-1.0)
	case "box-muller":
		u = rand(n); v = rand(n)
		x = np.sqrt(-2.0*np.log(u)) * np.cos(2.0*np.pi*v)
	# default random generator
	case _:
		x = rand(n)

# output samples, one per line, in scientific notation
for value in x:
	print("%.16e" % value)
