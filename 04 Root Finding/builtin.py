#!/usr/bin/env python
# demo of built-in root finders

#######################################################################

from math import *
from scipy.optimize import bisect, ridder, brentq, toms748, newton

#######################################################################

# find root of a function f(x) = 0
def f(x):
	return cos(x)

def df(x):
	return -sin(x)

# desired accuracy of argument value
epsilon = 1.0e-14

# bracketed interval
a = 0.0; fa = f(a)
b = 3.0; fb = f(b)

#######################################################################

c = bisect(f, a, b, xtol=epsilon); fc = f(c); print(c,fc)
c = ridder(f, a, b, xtol=epsilon); fc = f(c); print(c,fc)
c = brentq(f, a, b, xtol=epsilon); fc = f(c); print(c,fc)
c = toms748(f, a, b, xtol=epsilon); fc = f(c); print(c,fc)

c = newton(f, 0.5, df, tol=epsilon); fc = f(c); print(c,fc)
