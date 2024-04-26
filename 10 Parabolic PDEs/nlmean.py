#!/usr/bin/env python
# de-noise image using (fast) non-local means

#######################################################################

import numpy as np
from PIL import Image
from numpy.random import normal

#######################################################################

# high quality reference image (6x7 Kodak VPS III scan)
img = Image.open('vps.tif'); rgb = np.array(img)
src = 0.3*rgb[:,:,0] + 0.59*rgb[:,:,1] + 0.11*rgb[:,:,2]

# image degraded by additive white Gaussian noise
data = src + normal(scale=20.0, size=src.shape)
data = np.clip(data, 0.0, 255.0).astype(int)

#######################################################################

from numpy.linalg import norm
from scipy.ndimage import median_filter as median
from scipy.ndimage import gaussian_filter as gauss

# usual suspects (at about the same residual noise)
avg = gauss(data, 1); print(norm(avg-gauss(src, 1)))
med = median(data, 5); print(norm((med-median(src, 5))))

# image dimensions
h,w = src.shape
X,Y = np.meshgrid(np.arange(w)/(w-1), np.arange(h)/(h-1))

# A/B image split
def split(A, B):
	return np.where(X+Y<1, A, B)

#######################################################################

from joblib import Parallel, delayed

# feature space convolution parameters
sigma = 50.0; radius = 50.0; n = 256; jobs = 8

# luminance levels to sample
level = np.linspace(0.0, 255.0, n)

def nlmean(value):
	q = np.exp(-(data-value)**2/(2.0*sigma**2)); p = q*data
	return np.array([gauss(q, radius), gauss(p, radius)])

# convolution in feature space
stack = np.array(Parallel(n_jobs=jobs)(delayed(nlmean)(v) for v in level))

# filtered output placeholder
clean = np.zeros_like(data)

# this should really use interpolation with less slices
for i in range(h):
	for j in range(w):
		q,p = stack[data[i,j],:,i,j]
		clean[i,j] = p/q

#######################################################################

import matplotlib.pyplot as plt

fig = plt.figure(); fig.gca().set_aspect('equal')
plt.imshow(split(data,clean), vmin=0.0, vmax=255.0, cmap='gray', interpolation='none')
plt.contour(X+Y, levels=[1], colors='tab:red')

plt.show()
