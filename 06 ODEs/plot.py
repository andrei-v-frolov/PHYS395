#!/usr/bin/env python

# import libraries
import matplotlib

try:
	import astropy.io.fits as pyfits
except ImportError:
	import pyfits
import numpy as np

#import scipy.stats
from pylab import *
from mpl_toolkits.axes_grid1 import make_axes_locatable

# load trajectory data
hdu = pyfits.open('data.fit')

data = hdu[0].data
#green = hdu[1].data
#blue = hdu[2].data

ny,nx = data.shape

# data extent and meshes
x0 = hdu[0].header['CRVAL1']
y0 = hdu[0].header['CRVAL2']
dx = hdu[0].header['CDELT1']
dy = hdu[0].header['CDELT2']

x = x0 + arange(nx)*dx
y = y0 + arange(ny)*dy

X, Y = meshgrid(x, y)

hdu.close()

# create the figure
figure(figsize=(10,8), frameon=False)

# plot image data
c = matplotlib.colors.LinearSegmentedColormap.from_list("difference", ["blue", "white", "black"])
im = imshow(data, extent=(x[0],x[-1],y[0],y[-1]), origin='lower', aspect=dx/dy, cmap=c, interpolation='none')

# plot contours
#contour(X, Y, data, 32, cmap=cm.jet)

# make colorbar match the plot
divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", size="5%", pad=0.07)
plt.colorbar(im, cax=cax)

show()
