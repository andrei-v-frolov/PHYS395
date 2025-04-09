#!/usr/bin/env python
# solve Poisson equation in a grounded box using FFTs

#######################################################################

import numpy as np
from scipy.fft import dst, idst

#######################################################################

# number of grid points
nx = 256; ny = 256

# initialize uniform x and y grids
lx = 1.0; dx = 2.0*lx/nx
ly = 1.0; dy = 2.0*ly/ny

x = np.linspace(-lx+dx/2.0, lx-dx/2.0, nx)
y = np.linspace(-ly+dy/2.0, ly-dy/2.0, ny)

# wave numbers, note that DST-II SKIPS k=0
kx = np.pi/(2.0*lx) * np.arange(1,nx+1)
ky = np.pi/(2.0*ly) * np.arange(1,ny+1)

# 2D mesh LUTs
X,Y = np.meshgrid(x,y)
Kx,Ky = np.meshgrid(kx,ky)

#######################################################################

# smooth charge distribution centered at (x,y)
def charge(x=0, y=0, q=1, sigma=0.01):
	return q/(2.0*np.pi*sigma**2) * np.exp(-((X-x)**2+(Y-y)**2)/(2.0*sigma**2))

# solve Poisson equation L[phi] = -rho
def solve(rho):
	Rho = dst(dst(rho).T).T
	return idst(idst(Rho/(Kx*Kx+Ky*Ky)).T).T

# test function (and its Laplacian)
#f = np.exp(-(X*X+Y*Y)*4.5)
#df = 9.0*(9.0*(X*X+Y*Y)-2.0)*f
#residual = solve(df) + f

phi = solve(charge(0,0))
#phi = solve(charge(-1.0,-0.7,1.5,0.10) + charge(1.6,0.3,-2.0,0.15) + charge(-0.2,0.6,-2.0,0.20))

#######################################################################

import matplotlib.pyplot as plt

fig = plt.figure(); ax = fig.gca()

# warm color map
pot = plt.imshow(phi, cmap='OrRd', extent=[-lx,lx,-ly,ly], origin='lower', interpolation='none')
plt.colorbar()

levels = [np.max(phi)/2**(5-i) for i in range(5)]
style = {'levels': levels, 'cmap': 'Reds', 'norm': 'log', 'vmin': levels[0]/2}
iso = ax.contour(phi, extent=[-lx,lx,-ly,ly], **style)

'''
# symmetric color map
z0 = np.max(np.abs(phi))

pot = plt.imshow(phi, cmap='seismic', vmin=-z0, vmax=z0, extent=[-lx,lx,-ly,ly], origin='lower', interpolation='none')
plt.colorbar()

style = {'levels': 21, 'cmap': 'coolwarm'}
iso = ax.contour(phi, extent=[-lx,lx,-ly,ly], **style)
'''

#######################################################################

# interactive controls
from matplotlib.widgets import Slider

# make room for widgets
fig.subplots_adjust(left=0.2, bottom=0.15)

# user interface elements
x_slider = Slider(
    ax=fig.add_axes([0.21, 0.02, 0.55, 0.03]),
    valmin=-lx, valmax=lx, valinit=0.0,
    label='x'
)

y_slider = Slider(
    ax=fig.add_axes([0.07, 0.15, 0.0225, 0.73]),
    valmin=-ly, valmax=ly, valinit=0.0,
    orientation="vertical",
    label="y"
)

# function to be called anytime a slider value changes
def update(value):
	phi = solve(charge(x_slider.val,y_slider.val))
	pot.set_data(phi)
	global iso; iso.remove()
	iso = ax.contour(phi, extent=[-lx,lx,-ly,ly], **style)

# register update handler
x_slider.on_changed(update)
y_slider.on_changed(update)

plt.show()
