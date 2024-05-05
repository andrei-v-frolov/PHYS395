#!/usr/bin/env python
# synchrotron radiation demo for PHYS822

#######################################################################

from math import *
import numpy as np

#######################################################################

# number of field lines to plot
n = 32; pts = 128*3; dt = 0.03

# trajectory parametrization
phase = 0.0; beta = tanh(0.85)

# field line bundle
bundle = np.zeros([n,pts,4])

# charge trajectory (synchrotron)
def synchrotron(phase):
	c = cos(phase); s = sin(phase)
	return (s,c,beta*c,-beta*s)

# unit vector boosted with velocity v
def boosted(theta, vx, vy):
	v = sqrt(vx*vx+vy*vy)
	R = np.array([[vx,vy], [-vy,vx]])/v if v > 0 else np.eye(2)
	kx,ky = np.dot(R, [cos(theta),sin(theta)])
	return np.dot(R.T, [v+kx, sqrt(1.0-v*v)*ky])/(1.0+v*kx)

# march field line bundle forward
def march():
	global phase, bundle; phase += beta*dt
	x,y,vx,vy = synchrotron(phase)

	bundle = np.roll(bundle, 1, axis=1)
	bundle[:,:,0:2] += bundle[:,:,2:4]*dt

	for i in range(n):
		kx,ky = boosted(2.0*np.pi*i/n, vx, vy)
		bundle[i,0] = [x,y,kx,ky]

# initial transient
for i in range(pts): march()

#######################################################################

import matplotlib.pyplot as plt

fig = plt.figure(); ax = fig.gca()

# circle traced by the charge
ax.add_patch(plt.Circle((0, 0), 1.0, color='tab:orange', linewidth=3, fill=False))

# field lines and charge location
field = [plt.plot(bundle[i,:,0], bundle[i,:,1], color='tab:blue')[0] for i in range(n)]
charge, = plt.plot([0], [1], "ro")

ax.set_aspect('equal')
plt.xlim([-7,7])
plt.ylim([-5,5])

#######################################################################

# interactive controls
from matplotlib.widgets import Slider, RadioButtons

# make room for widgets
fig.subplots_adjust(bottom=0.10)

psi_slider = Slider(
    ax=fig.add_axes([0.125, 0.02, 0.775, 0.03]),
    valmin=0.0, valmax=1.5, valinit=0.75,
    label='ψ'
)

# function to be called anytime a slider value changes
def update(value):
	global beta; beta = tanh(psi_slider.val)

# register update handler
psi_slider.set_val(atanh(beta))
psi_slider.on_changed(update)

#######################################################################

import matplotlib.animation as animation

# called to advance animation to next frame
def animate(i):
	march(); x,y,vx,vy = synchrotron(phase); charge.set_data([x],[y])
	for i in range(n): field[i].set_data(bundle[i,:,0], bundle[i,:,1])

animation = animation.FuncAnimation(fig, animate, frames=1000, interval=1000.0/60)
#animation.save('synchrotron.mp4')
plt.show()
