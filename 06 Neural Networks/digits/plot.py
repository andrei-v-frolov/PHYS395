#!/usr/bin/env python
# render training and test MNIST handwritten data images

# parse arguments
from sys import argv
file = argv[1] if len(argv) > 1 else None

# import libraries
import idx
import matplotlib
matplotlib.use('macosx' if file is None else 'PDF')

from pylab import *

# load mnist test data
test = idx.open('data/t10k-images-idx3-ubyte.gz')
image = idx.grid(test, 100, 100)
#image = idx.grid(test, 125, 80)
#image = idx.grid(test, 200, 50)

#test = idx.open('data/train-images-idx3-ubyte.gz')
#image = idx.grid(test, 300, 200)

# create the figure
figure(figsize=(10,8), frameon=False); gradient = ["white", "black"]
c = matplotlib.colors.LinearSegmentedColormap.from_list("difference", gradient)
imshow(image, origin='upper', cmap=c, vmin=0.0, vmax=255.0, interpolation='none')

tick_params(left=False, right=False, labelleft=False, labelbottom = False, bottom=False)

if not(file is None): plt.savefig(file, bbox_inches='tight', pad_inches=0.0, transparent=True)
show()