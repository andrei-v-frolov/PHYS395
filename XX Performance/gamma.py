#!/usr/bin/env python
# effective gamma of matplotlib colormaps

#######################################################################

import numpy as np

#######################################################################

# test target resolution
n = 1024; block = 8

# test target patterns
X,Y = np.meshgrid(range(n),range(n))

checker = (X//block + Y//block) % 2
cross = (X//(n>>1) + Y//(n>>1)) % 2
disk = (X-n/2)**2 + (Y-n/2)**2 < 1.5*(n/4)**2
solid = np.logical_xor(cross,disk)

# test target alternating between solid and stippled colors
def target(gamma=1.0):
	return np.where(solid, 0.5**gamma, checker)

#######################################################################

import matplotlib.pyplot as plt

fig = plt.figure(); ax = fig.gca()

# test target render
pattern = plt.imshow(checker, cmap='Greys', origin='lower', interpolation='none')
plt.colorbar()

ax.set_aspect('equal')

#######################################################################

# interactive controls
from matplotlib.widgets import Slider

# make room for widgets
fig.subplots_adjust(bottom=0.13)

# user interface elements
gamma_slider = Slider(
    ax=fig.add_axes([0.18, 0.02, 0.57, 0.03]),
    valmin=0.5, valmax=2.5, valinit=1.5,
    label='ɣ'
)

# function to be called anytime a slider value changes
def update(value):
	pattern.set_data(target(gamma_slider.val))

# register update handler
gamma_slider.on_changed(update)

update(None)
plt.show()
