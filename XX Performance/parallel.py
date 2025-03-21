#!/usr/bin/env python
# parallel versus vectorized performance demo

#######################################################################

from time import perf_counter as now
from joblib import Parallel, delayed
from joblib import parallel_config

#######################################################################

from math import *
import numpy as np

# number of items to process
n = 100000

#######################################################################
# vectorize light-weight operation using NumPy
#######################################################################

t1 = now()
squares = np.fromiter((i**2 for i in range(n)), dtype='float')
t2 = now(); print("NumPy array fromiter(): ", t2-t1)

t1 = now()
squares = np.sqrt(np.arange(n)**2)
t2 = now(); print("NumPy native arange():  ", t2-t1)

t1 = now()
np.vectorize(sqrt)(squares)
t2 = now(); print("NumPy vectorize(sqrt): ", t2-t1)

t1 = now()
np.sqrt(squares)
t2 = now(); print("NumPy np.sqrt():       ", t2-t1)

#######################################################################
# execute in parallel threads using Joblib
#######################################################################

with parallel_config(backend='threading', n_jobs=16):
	# worker pool creation
	t1 = now()
	parallel = Parallel()
	t2 = now(); print("Joblib worker pool:    ", t2-t1)

	# time first invocation
	t1 = now()
	parallel(delayed(sqrt)(i**2) for i in range(n))
	t2 = now(); print("Joblib first batch:    ", t2-t1)

	# time second invocation
	t1 = now()
	parallel(delayed(sqrt)(i**2) for i in range(n))
	t2 = now(); print("Joblib second batch:   ", t2-t1)
