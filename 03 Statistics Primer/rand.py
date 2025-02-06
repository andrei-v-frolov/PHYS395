#!/usr/bin/env python
# random number generator demo

#######################################################################

# using system random number generator
from random import random

# draw random numbers one by one
for i in range(100):
	print(random())

#######################################################################

# bulk random number generator in NumPy
import numpy as np
from numpy.random import seed, rand

# number of samples to draw
n = 1024

# one can seed the generator for repeatable runs
# seed(1)

# generate random samples
x = rand(n)

# simple mean and variance estimators
mean = np.sum(x)/n
variance = np.sum((x-mean)**2)/(n-1)

# note that estimators are random variables!
print(mean, variance)
