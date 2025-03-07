#!/usr/bin/env python
# MCMC sampler of fit parameters for Assignment #4

#######################################################################

import numpy as np
from scipy.special import loggamma
from numpy.random import normal, uniform

#######################################################################

# load data to fit
data = np.loadtxt('periodic.dat')
x = data[:,0]; y = data[:,1]

# model to be fitted
def f(c,x):
	k = 2.0*np.pi*np.arange(len(c)-1)
	return np.exp(c[1:] @ np.cos(np.outer(k,x))) + c[0]

'''
# logarithm of chi^2 likelihood tail
def loglikelihood(c):
	k,_ = data.shape; k -= len(c)
	z = np.sum((y-f(c,x))**2)/k
	return (k/2)*(1.0-z+np.log(z))
'''

# logarithm of chi^2 PDF (up to a constant)
def loglikelihood(c):
	k,_ = data.shape; k -= len(c)
	z = np.sum((y-f(c,x))**2)/2
	return (k/2-1)*np.log(z)-z-loggamma(k/2)

#######################################################################

# number of coefficients and initial guess
n = 4; c = np.zeros(n); sigma = 0.03

# MCMC step using Metropolis-Hastings algorithm
def step(x):
	y = x + normal(scale=sigma, size=n)
	delta = loglikelihood(y)-loglikelihood(x)
	return np.where(uniform() > np.exp(delta), x, y)

#for i in range(1000): c = step(c)

#######################################################################

import matplotlib.pyplot as plt

fig = plt.figure()

plt.plot(x, y, "+")
fit, = plt.plot(x, f(c,x), "r-")
plt.xlim([0,1])

#######################################################################

import matplotlib.animation as animation

# called to advance animation to next frame
def animate(i):
	global c; c = step(c)
	fit.set_data(x, f(c,x))

animation = animation.FuncAnimation(fig, animate, frames=2000, interval=1000/60)
#animation.save('mcmc.mp4')
plt.show()
