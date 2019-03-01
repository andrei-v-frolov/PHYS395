#!/usr/bin/env python

import matplotlib
matplotlib.use('macosx')

import scipy.integrate
from pylab import *

def f(y,t):
	return y[1],-y[0]**3

t = arange(0.0,100.0,1.0) * 7.416298709205487673735401388781040185
y = scipy.integrate.odeint(f, [1.0,0.0], t, atol=1.0e-12)
e = y[:,0]**4/4.0 + y[:,1]**2/2.0 - 0.25

figure(figsize=(10,3), frameon=False)

plot(t, e, 'r-')

xlabel('t')
ylabel('y')

show()

#savefig("smpout.pdf", format="pdf")
