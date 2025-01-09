#!/usr/bin/env python
# Fibonacci numbers using a list

# initialize the sequence
f = [0,1]

for i in range(2,300):
	f.append(f[i-2]+f[i-1])

print(f)