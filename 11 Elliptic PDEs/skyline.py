#!/usr/bin/env python
# potential between clouds and Vancouver skyline

#######################################################################

import numpy as np
from PIL import Image
from multigrid import grid

#######################################################################

# load skyline mask
img = Image.open('mask.png'); rgb = np.array(img)
mask = (0.3*rgb[:,:,0] + 0.59*rgb[:,:,1] + 0.11*rgb[:,:,2]) < 220

# initialize boundary conditions
h,w = mask.shape
seed = np.ones([h,w])
seed[0:h>>1,:] = -1.0

phi = np.where(mask, seed, 0.0)

# find the potential
grid(phi,1.0,mask=mask,bc='reflect').solve()

# compute the electric field strength
gy,gx = np.gradient(phi,1.0)
E = np.sqrt(gx*gx+gy*gy)

#######################################################################

import matplotlib.pyplot as plt
from matplotlib import colors

overlay = colors.LinearSegmentedColormap.from_list("mask", ['#00000000', 'black'])

fig = plt.figure(); ax = fig.gca()
#plt.imshow(E, cmap='seismic', interpolation='none')
plt.imshow(E, cmap='twilight', norm='log', vmin=1.0e-3, interpolation='none')
plt.colorbar(fraction=0.0347, pad=0.05)
#plt.imshow(mask, cmap=overlay, interpolation='none')

plt.show()
