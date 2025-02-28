#!/usr/bin/env python
# demo of built-in SciPy optimization methods

#######################################################################

from scipy.optimize import minimize

#######################################################################

# tilt parameter
mu = 0.1

# tilted Mexican hat (testing convergence in winding valley)
def f(u):
	x,y = u; return (x*x+y*y - 1.0)**2/4.0 - mu*x

#######################################################################

# starting point
start = [0.0,2.0]

# if not supplied, gradient is evaluated by finite difference
result = minimize(f, start); print("SciPy default", result)
result = minimize(f, start, method='BFGS'); print("BFGS", result)
result = minimize(f, start, method='L-BFGS-B'); print("L-BFGS-B", result)
result = minimize(f, start, method='CG'); print("Conjugate gradient", result)
result = minimize(f, start, method='Powell'); print("Powell's method", result)
result = minimize(f, start, method='Nelder-Mead'); print("Nelder-Mead method", result)

