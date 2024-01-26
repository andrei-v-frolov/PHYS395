#!/usr/bin/env python
# compute conditioning number for sampled orthogonal polynomial basis

#######################################################################

import numpy as np
from numpy.polynomial.legendre import legvander, legval
from numpy.polynomial.chebyshev import chebvander, chebval
from numpy.linalg import solve, eig, svd

#######################################################################

# number of coefficients
n = 30

# evaluation grid - uniform
x = np.linspace(-1.0, 1.0, n)

# evaluation grid - Chebyshev
#theta = np.linspace(np.pi, 0.0, n)
#x = np.cos(theta)

# function to approximate
f = np.exp(-x*x*4.5)

# coefficients satisfy A*c = f
A = legvander(x,n-1)
S = svd(A, compute_uv=False)

# compute conditioning number
print(S[0]/S[-1])

# if S was *not* sorted, we could do
#S.sort(kind='stable'); print("%.2E" % (S[-1]/S[0]))
