#!/usr/bin/env python
# gradient descent in 1D

#######################################################################

from math import *
import numpy as np

#######################################################################

# function to be minimized
def f(x):
	return x**4/4.0

def df(x):
	return x**3

def ddf(x):
	return 3.0*x*x

# number of iterations
n = 200

# step size (aka learning rate)
eta = 0.1

# initial guess
x = 2.0; g = df(x)

# optimization history
history = np.zeros(n)

# fixed number of iterations
for i in range(0,n):
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
	print(x, f(x))

#######################################################################

import matplotlib.pyplot as plt

plt.plot(history,'r-')
plt.yscale('log')

plt.show()
