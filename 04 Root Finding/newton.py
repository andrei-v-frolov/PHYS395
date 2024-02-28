#!/usr/bin/env python
# find a root of a function using Newton's method

#######################################################################

from math import *

#######################################################################

# find root of a function f(x) = 0
def f(x):
	return cos(x)

def df(x):
	return -sin(x)

# initial guess
x = 0.5

# Newton's method update
for i in range(0,8):
	x -= f(x)/df(x)
	print(x,f(x))
