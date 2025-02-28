#!/usr/bin/env python
# find a minimum by gradient descent in 1D

#######################################################################

from math import *
import numpy as np

#######################################################################

# function to be minimized
def f(x):
	return 1.0-cos(x)

def df(x):
	return sin(x)

def ddf(x):
	return cos(x)

#######################################################################

# number of iterations
n = 200

# step size (aka learning rate)
eta = 0.1

# initial guess
x = 1.0; g = df(x)

# optimization history
history = np.zeros(n)

# fixed number of iterations
for i in range(n):
	# fixed step
	#x -= copysign(eta,df(x))
	
	# gradient descent
	#x -= eta*df(x)
	
	# Newton's step
	#x -= df(x)/ddf(x)
	
	# Barzilai-Borwein
	dx = -eta*g; x += dx
	dg = df(x)-g; g += dg
	eta = dx/dg
	
	# log history
	history[i] = x
	print(i, x, f(x))

#######################################################################

import matplotlib.pyplot as plt

plt.plot(history,'r-')
plt.yscale('log')

plt.show()
