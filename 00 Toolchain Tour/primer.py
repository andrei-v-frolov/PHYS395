#!/usr/bin/env python
# Python language features primer

#######################################################################
# intrinsic data types (note that Python is *not* strongly typed)
#######################################################################

# integers (64-bit natively, but arbitrary precision is supported!)
i = 1; j = i<<72; print(i, j)

# floats (double precision, 64 bits, about 16 significant digits)
f = 1.0; print(f, ('%+.16f' % f))

# floats are fixed precision and *do* overflow!
g = 1.0e300*1.0e300; print(g)

# complex numbers are natively supported
z = 0.0 + 1.0j; print(z*z)

# string literals have many forms, simplest is
s = "Hello, World!"; print(s)

# list is an ordered set of values (of the same type)
l = [1, 2, 3]; print(l, l[0], l[1], l[-1])

# dictionary is unordered lookup table mapping keys to values
color = {'apple': 'red', 'orange': 'orange', 'banana': 'yellow'}
print("apple is " + color['apple'])

#######################################################################
# loading modules and doing simple calculations
#######################################################################

from math import *

print(sin(1.0), exp(1.0))

#######################################################################
# NumPy arrays (multidimensional, support slicing, Fortran syntax)
#######################################################################

import numpy as np

# simple 1D array from a list
a = np.array([1, 2, 3, 4, 5]); print(a)

# array slicing
print(a[2:], a[:-2], a[::2])

# multiple constructors are available
a = np.zeros([2,3]); print(a)
a = np.identity(3); print(a)
a = np.diag([1,2,3]); print(a)
a = np.eye(2,3); print(a)

# array slicing
print(a[0,:], a[:,1])

#######################################################################
# flow control: conditionals, ternary operator, loops and iterators
#######################################################################

# conditional execution
x = 15

if x < 3:
	print(x, "is a few")
elif x < 10:
	print(x, "is a handful")
else:
	print(x, "is a lot")

# conditional expressions
big = True
size = 16.0 if big else 8.0
print(size)

# iterating over collections
for i in range(15):
	print(i)

for i in l:
	print(i)

for k,v in color.items():
	print(k,v)

# conditional loop
i = 0
while i < 15:
	i += 1; print(i)

#######################################################################
# functions and closures
#######################################################################

# named arguments and default values
def scale(x, mu=0.0, sigma=1.0):
	return asinh((x-mu)/sigma)

print(scale(1.0), scale(1.0, sigma=3.0))

# closures are anonymous functions
list = [13,2,17,21,3,0,15]
print(sorted(list, key=lambda x: +x))
print(sorted(list, key=lambda x: -x))
