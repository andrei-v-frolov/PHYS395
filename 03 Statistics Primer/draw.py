#!/usr/bin/env python
# draw a specified number of IID samples from a distribution
# run as: python draw.py [n:samples] [seed:value]

#######################################################################

# access command line arguments
from sys import argv

# default argument values
n = 1024; init = None

# parse command line arguments (if supplied)
for arg in argv[1:]:
	match tuple(arg.split(":")):
		case ("n",s):
			n = int(s)
		case ("seed",s):
			init = int(s)
		case _:
			...

#######################################################################

# bulk random number generators from NumPy
import numpy as np
from numpy.random import seed, rand

# seed random sequence
seed(init)

# generate random samples
x = rand(n)

# output samples, one per line, in scientific notation
for value in x:
	print("%.16e" % value)
