#!/usr/bin/env python
# interactive demo of ballistic motion

#######################################################################

from math import *
import numpy as np
from scipy.integrate import solve_ivp

#######################################################################

# simulation parameters
g = 9.81; alpha = 0.0; beta = 0.0

# initial conditions
v = 10.0; angle = 45.0

# equations of motion for projectile with air drag
def f(t,state):
	x,y,vx,vy = state
	drag = alpha + beta*sqrt(vx*vx+vy*vy)
	return [vx,vy,-drag*vx,-g-drag*vy]

# stopping condition (projectile hits the ground)
def stop(t,state):
	x,y,vx,vy = state
	return y

# stopping condition attributes (terminate on impact)
stop.terminal = True
stop.direction = -1.0

#######################################################################

# integrate projectile trajectory
def trajectory(v, angle, n=256):
	phi = (pi/180.0) * angle
	t = np.linspace(0.0, 2.0*v/g, n)
	ic = [0.0,0.0,v*cos(phi),v*sin(phi)]
	soln = solve_ivp(f, [t[0],t[-1]], ic, t_eval=t, events=stop)
	return (soln.y[0],soln.y[1])

x,y = trajectory(v,angle)

#######################################################################

import matplotlib.pyplot as plt

# plot the trajectory
fig,main = plt.subplots()
path, = main.plot(x,y,'r-')

# set maximal bounds
h = v*v/(2.0*g)
main.set_xlim([0,2.04*h])
main.set_ylim([0,1.02*h])

#######################################################################

# interactive controls
from matplotlib.widgets import Slider

# keep only n significant digits
def digits(x, n=1):
	return round(x, n-int(ceil(log10(abs(x)))))

# make room for widgets
fig.subplots_adjust(left=0.2, bottom=0.22)

# user interface elements
v_slider = Slider(
    ax=fig.add_axes([0.07, 0.22, 0.0225, 0.66]),
    valmin=0, valmax=v, valinit=v,
    orientation="vertical",
    label="v [m/s]"
)

angle_slider = Slider(
    ax=fig.add_axes([0.2, 0.12, 0.7, 0.03]),
    valmin=0.0, valmax=90.0, valinit=angle,
    label='angle [deg] '
)

alpha_slider = Slider(
    ax=fig.add_axes([0.2, 0.07, 0.7, 0.03]),
    valmin=0.0, valmax=digits(3.0*g/v), valinit=alpha,
    label='alpha [1/s] '
)

beta_slider = Slider(
    ax=fig.add_axes([0.2, 0.02, 0.7, 0.03]),
    valmin=0.0, valmax=digits(10.0*g/(v*v)), valinit=beta,
    label='beta [1/m] '
)

# function to be called anytime a slider value changes
def update(value):
	global alpha, beta
	alpha = alpha_slider.val; beta = beta_slider.val
	x,y = trajectory(v_slider.val, angle_slider.val)
	path.set_data(x,y); fig.canvas.draw_idle()

# register update handler
v_slider.on_changed(update)
angle_slider.on_changed(update)
alpha_slider.on_changed(update)
beta_slider.on_changed(update)

plt.show()
