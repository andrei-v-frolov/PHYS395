#!/usr/bin/env python
# using high-resolution timer for profiling

#######################################################################

from time import perf_counter as now

#######################################################################

# initialize accumulator
a = 0

t1 = now()

# block of code to be timed
for i in range(0,1024*100): a += i

t2 = now()

print(t2-t1)