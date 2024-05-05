#!/usr/bin/env python
# local entropy estimator for an image

#######################################################################

import numpy as np
from PIL import Image
from scipy.ndimage import median_filter as median
from skimage.filters.rank import entropy
from skimage.morphology import disk

#######################################################################

# estimation radius
n = 11

# realistic AI-augmented landscape image
img = Image.open('trecime.jpg'); rgb = np.array(img)
src = 0.3*rgb[:,:,0] + 0.59*rgb[:,:,1] + 0.11*rgb[:,:,2]

# extract noise component
src -= median(src,n)

# local entropy estimator
out = entropy((128+3*src).astype(np.uint8), disk(n))

# downsample source and output images
h,w = src.shape; img = np.array(img.resize((w//n,h//n), resample=Image.LANCZOS))
out = np.array(Image.fromarray(out).resize((w//n,h//n), resample=Image.LANCZOS))

#######################################################################

# scale entropy and output geometry
h,w = out.shape; out *= 255.0/np.max(out)
X,Y = np.meshgrid(np.arange(w)/(w-1), np.arange(h)/(h-1))

# A/B image split
for i in range(3):
	img[:,:,i] = np.where(X+Y<1, img[:,:,i], out)

#######################################################################

import matplotlib.pyplot as plt

fig = plt.figure(); fig.gca().set_aspect('equal')
#plt.imshow(src, vmin=0.0, vmax=255.0, cmap='gray', interpolation='none')
plt.imshow(img, interpolation='none')
plt.contour(X+Y, levels=[1], colors='tab:red')

plt.show()
