#!/usr/bin/env python
# render success rate on test MNIST handwritten data images

#######################################################################

# parse arguments
from sys import argv
file = argv[1] if len(argv) > 1 else None

#######################################################################

# import PyTorch libraries
import torch
from model import *

# current device for testing
print(f"Using {device} device")

# load trained model
model = NeuralNetwork().to(device)
model.load_state_dict(torch.load("model.pth"))
model.eval()
print(model)

#######################################################################

from dataset import *

# test data and success bitmap
data = HandwrittenDigitsDataset('t10k')

with torch.no_grad():
	pred = model(data.image.to(device))
	success = (pred.argmax(-1).cpu() == data.label)

# tile the data for plotting
nx = 125; ny = 80
image = idx.tile(data.image, nx, ny)
success = idx.tile(success.reshape((-1,1,1)), nx, ny)

#######################################################################

# import matplotlib libraries
import matplotlib, matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap as colormap

# if saving to file, use PDF backend
if file is not None:
	matplotlib.use('PDF')

# create the figure
plt.figure(figsize=(10,8), frameon=False)

# classification success image
s = colormap.from_list("success", ["red", "white"])
plt.imshow(success, origin='upper', extent=[0,nx,0,ny], cmap=s, vmin=0.0, vmax=1.0, interpolation='none')

# test data overlay image
d = colormap.from_list("overlay", ["#00000000", "black"])
plt.imshow(image, origin='upper', extent=[0,nx,0,ny], cmap=d, vmin=0.0, vmax=1.0, interpolation='none')

# no ticks and tight framing
plt.tick_params(left=False, right=False, labelleft=False, labelbottom = False, bottom=False)
plt.tight_layout()

# show or save the figure
plt.show() if file is None else plt.savefig(file, bbox_inches='tight', pad_inches=0.0, transparent=True)
