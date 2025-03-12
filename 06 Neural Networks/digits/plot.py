#!/usr/bin/env python
# render training and test MNIST handwritten data images

#######################################################################

# parse arguments
from sys import argv
file = argv[1] if len(argv) > 1 else None

#######################################################################

import idx

# load mnist test data
test = idx.open('data/t10k-images-idx3-ubyte.gz')
image = idx.tile(test, 125, 80)

#train = idx.open('data/train-images-idx3-ubyte.gz')
#image = idx.tile(train, 300, 200)

#######################################################################

# import matplotlib libraries
import matplotlib, matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap as colormap

# if saving to file, use PDF backend
if file is not None:
	matplotlib.use('PDF')

# create the figure
plt.figure(figsize=(10,8), frameon=False)

# bitmap image with custom colormap
cmap = colormap.from_list("B&W", ["white", "black"])
plt.imshow(image, origin='upper', cmap=cmap, vmin=0.0, vmax=255.0, interpolation='none')

# no ticks and tight framing
plt.tick_params(left=False, right=False, labelleft=False, labelbottom = False, bottom=False)
plt.tight_layout()

# show or save the figure
plt.show() if file is None else plt.savefig(file, bbox_inches='tight', pad_inches=0.0, transparent=True)
