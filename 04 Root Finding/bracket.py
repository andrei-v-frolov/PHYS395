#!/usr/bin/env python
# find a bracketed root by bisection (and friends)

#######################################################################

from math import *

#######################################################################

# find root of a function f(x) = 0
def f(x):
	return cos(x)

# desired accuracy of argument and function values
epsilon = 1.0e-14; delta = 1.0e-16

# bracketed interval
a = 0.0; fa = f(a)
b = 2.0; fb = f(b)

# check that it is indeed bracketed
assert fa*fb < 0, f'Interval [{a},{b}] does not bracket a root of f(x)!'

# shrink the interval to specified size
while fabs(b-a) > epsilon:
	# straight bisection
	c = (a+b)/2.0; fc = f(c)
	
	# Ridder's method
	c += (c-a)*copysign(1.0,fa-fb)*fc/sqrt(fc*fc-fa*fb); fc = f(c)
	
	# false position method
	#c = (a*fb-b*fa)/(fb-fa); fc = f(c)
	
	# bail out if we hit the root or get stuck
	if (abs(fc) < epsilon or c == a or c == b): break
	
	# output progress report
	print(c, fc)
	
	# update bracketed interval
	if fa*fc < 0.0: b,fb = c,fc
	if fc*fb < 0.0: a,fa = c,fc

'''
# output final value
print(c, fc)
'''
