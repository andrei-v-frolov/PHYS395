#!/usr/bin/env python
# compute conditioning number for least square fit

#######################################################################

import numpy as np
from numpy.polynomial.legendre import legvander, legval
from numpy.polynomial.chebyshev import chebvander, chebval
from numpy.linalg import solve, eig, svd

#######################################################################

# number of coefficients
n = 10; pts = 10

# evaluation grid - uniform
x = np.linspace(-1.0, 1.0, pts)

# evaluation grid - Chebyshev
#theta = np.linspace(np.pi, 0.0, n)
#x = np.cos(theta)

# coefficients satisfy B*c = y
B = legvander(x,n-1)
A = np.matmul(B.T,B)
S = svd(A, compute_uv=False)

# compute conditioning number
print(S[0]/S[-1])

# if S was *not* sorted, we could do
#S.sort(kind='stable'); print("%.2E" % (S[-1]/S[0]))
