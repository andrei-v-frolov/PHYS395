#!/usr/bin/env python
# Fibonacci numbers using ring buffer

# initialize the sequence
a = 0
b = 1

for i in range(2,100):
	c = a+b; print(c)
	a = b+c; print(a)
	b = c+a; print(b)
