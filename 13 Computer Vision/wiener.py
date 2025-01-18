#!/usr/bin/env python
# de-noise image using simple Wiener filter (aka H^1 regularization)

#######################################################################

import numpy as np
from PIL import Image
from scipy.fft import dct, idct
from numpy.random import normal

#######################################################################

# high quality reference image (6x7 Kodak VPS III scan)
img = Image.open('vps.tif'); rgb = np.array(img)
src = 0.3*rgb[:,:,0] + 0.59*rgb[:,:,1] + 0.11*rgb[:,:,2]

# image degraded by additive white Gaussian noise
data = src + normal(scale=20.0, size=src.shape)
data = np.clip(data, 0.0, 255.0).astype(int)

#data = np.array(Image.open('terminal.gif'))

#######################################################################

# number of grid points
ny,nx = data.shape

# wave numbers, note that DCT-II includes k=0
kx = np.pi/(2.0*nx) * np.arange(0,nx)
ky = np.pi/(2.0*ny) * np.arange(0,ny)

# 2D mesh LUTs
Kx,Ky = np.meshgrid(kx,ky)
K2 = Kx*Kx+Ky*Ky

#######################################################################

# solve massive Poisson equation L[phi] - mu^2*phi = -rho
def solve(rho, mu=1.0):
	Rho = dct(dct(rho).T).T
	Phi = Rho/(K2 + mu*mu)
	if not np.isfinite(Phi[0,0]): Phi[0,0] = 0.0
	return idct(idct(Phi).T).T

phi = solve(data)

#######################################################################

import matplotlib.pyplot as plt

fig = plt.figure(); ax = fig.gca()

out = plt.imshow(phi, origin='upper', cmap='gray', vmin=0.0, vmax=255.0, interpolation='none')
plt.colorbar()

#######################################################################

# interactive controls
from matplotlib.widgets import Slider

# make room for widgets
fig.subplots_adjust(bottom=0.15)

# user interface elements
scale = Slider(
    ax=fig.add_axes([0.21, 0.02, 0.55, 0.03]),
    valmin=-1.0, valmax=1.0, valinit=0.0,
    label='log₁₀ρ'
)

# function to be called anytime a slider value changes
def update(value):
	mu = 10.0**(-scale.val)
	phi = solve(mu*mu*data, mu=mu)
	out.set_data(phi)

# register update handler
scale.on_changed(update)

plt.show()
