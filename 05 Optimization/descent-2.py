#!/usr/bin/env python
# gradient descent in 2D

#######################################################################

import numpy as np

#######################################################################

# tilt parameter
mu = 0.1

# tilted Mexican hat (testing convergence in winding valley)
def f(x,y):
	return (x*x+y*y - 1.0)**2/4.0 - mu*x

def df(x,y):
	return np.array([(x*x+y*y - 1.0)*x - mu, (x*x+y*y - 1.0)*y])

def ddf(x,y):
	return np.array([[3.0*x*x+y*y - 1.0, 2.0*x*y], [2.0*x*y, 3.0*x*x+y*y - 1.0]])

#######################################################################

# number of iterations
n = 200

# step size (aka learning rate)
eta = 0.1

# initial guess
x = [0.0,2.0]; g = df(*x)

# optimization history
history = np.zeros([2,n])

# fixed number of iterations
for i in range(0,n):
	# fixed step
	#x -= eta*df(*x)/np.linalg.norm(df(*x))
	
	# gradient descent
	#x -= eta*df(*x)
	
	# Newton's step
	#x -= np.linalg.solve(ddf(*x), df(*x))
	
	# Barzilai-Borwein
	dx = -eta*g; x += dx
	dg = df(*x)-g; g += dg
	eta = min(np.dot(dx,dx)/np.dot(dx,dg), np.dot(dx,dg)/np.dot(dg,dg))
	
	# log history
	history[:,i] = x
	print(x, f(*x))

#######################################################################

import matplotlib.pyplot as plt

# domain to be plotted
l = 1.5; pts = 1024

# function being minimized
grid = np.linspace(-l,l,pts)
X,Y = np.meshgrid(grid,grid)
F = np.vectorize(f)(X,Y)

# scale function for plotting
F = F/(mu*mu+F*F)**0.375

# plot optimization history
plt.plot(history[0], history[1], "r.-")
plt.imshow(F, extent=[-l,l,-l,l], cmap='twilight')
#plt.contourf(F, extent=[-l,l,-l,l], levels=30, cmap='twilight')
plt.colorbar()

# make sure aspect is 1:1
plt.gca().set_aspect('equal')

plt.show()