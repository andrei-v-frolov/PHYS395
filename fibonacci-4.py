#!/usr/bin/env python
# Fibonacci numbers using tuples

# initialize the sequence
a = 0
b = 1

for i in range(2,100):
	a,b = b,a+b
	print(b)
