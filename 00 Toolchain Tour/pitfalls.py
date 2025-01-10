#!/usr/bin/env python
# common pitfalls in Python numerical code

#######################################################################

import numpy as np

#######################################################################

# floating point vs integer division (behavior differs in v2 and v3!!!)
print(4/3, 4//3)

# max integer size in intrinsic vs NumPy types
i = 1 << 72; n = np.array([1]) << 72; print(i,n[0])

# limited float precision can lead to unexpected results
n = 1000000
a = np.array([0], dtype='int64')
b = np.array([0], dtype='float32')
for i in range(n):
	a += n-i; b += n-i
print(a[0],b[0])

# exception handling by NumPy is different!
print(np.zeros(1)/0.0)
print(0.0/0.0)
