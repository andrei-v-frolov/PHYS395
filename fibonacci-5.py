#!/usr/bin/env python
# Fibonacci numbers using recursion

def fibonicci(a,b):
	c = a+b; print(c)
	fibonicci(b,c)

fibonicci(0,1)