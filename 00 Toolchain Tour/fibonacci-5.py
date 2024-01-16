#!/usr/bin/env python
# Fibonacci numbers using recursion

def fibonacci(a,b):
	c = a+b; print(c)
	fibonacci(b,c)

fibonacci(0,1)