#!/usr/bin/env python

# import libraries
import matplotlib

try:
	import astropy.io.fits as pyfits
except ImportError:
	import pyfits
import numpy as np

from pylab import *
from scipy.signal import medfilt
from mpl_toolkits.axes_grid1 import make_axes_locatable

# load trajectory data
hdu = pyfits.open('data.fit')

data = hdu[0].data

ny,nx = data.shape

# antialias using median filter
osy = 1; osx = 1

if (osy != 1 or osx != 1):
	aakernel = [(osy & (~1)) + 1, (osx & (~1)) + 1]
	data = medfilt(data, aakernel)[::osy,::osx]; ny,nx = data.shape
	print("Antialiased with median kernel %s; output size is (%i,%i) pixels" % (aakernel,nx,ny))

# data extent and meshes
x0 = hdu[0].header['CRVAL1']
y0 = hdu[0].header['CRVAL2']
dx = hdu[0].header['CDELT1']
dy = hdu[0].header['CDELT2']

x = x0 + arange(nx)*dx*osx
y = y0 + arange(ny)*dy*osy

X, Y = meshgrid(x, y)

hdu.close()

# create the figure
figure(figsize=(10,8), frameon=False)

# plot image data
c = matplotlib.colors.LinearSegmentedColormap.from_list("difference", ["blue", "white", "black"])
im = imshow(data, extent=(x[0],x[-1],y[0],y[-1]), origin='lower', aspect=(osx*dx)/(osy*dy), cmap=c, interpolation='none')

# plot contours
#contour(X, Y, data, 32, cmap=cm.jet)

# make colorbar match the plot
divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", size="5%", pad=0.07)
plt.colorbar(im, cax=cax)

show()
