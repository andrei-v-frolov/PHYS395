#!/usr/bin/env python
# visualize neural network operation and link weights

#######################################################################

# parse arguments
from sys import argv
file = argv[1] if len(argv) > 1 else None

#######################################################################

from grid import *

# geometry of network layers for plotting
layer = [
	grid((28,28), dx=(1,-1), s=18),
	grid((11,11), dx=(2.5,-2.5), dy=(0,2.5), center=(36,0), s=120),
	grid((5,5), dx=(3.5,-1), dy=(0,3.5), center=(65,0), s=250),
	grid((10,1), dy=(0,5.5), center=(84,0), s=800)
]

#######################################################################

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection as lines
from matplotlib.colors import CenteredNorm, LinearSegmentedColormap as colormap

# if saving to file, use PDF backend
if file is not None:
	matplotlib.use('PDF')

# create the figure (1080p, prevent scaling)
fig = plt.figure(figsize=(48/5,27/5), frameon=False)
ax = fig.gca(); ax.set_aspect('equal')

# plot neural network nodes
cmap = colormap.from_list('neuron', ['lightsteelblue', 'black'])
nodes = [plt.scatter(*l.array().T, c=np.zeros(l.count), cmap=cmap, vmin=0.0, vmax=1.0, zorder=2, **l.kwargs) for l in layer]

# output node labels
for i,x in enumerate(layer[-1]):
	plt.text(x[0], x[1]-0.2, i, color='white', weight='bold', fontsize=18, horizontalalignment='center', verticalalignment='center')

# no ticks and tight framing
plt.tick_params(left=False, right=False, labelleft=False, labelbottom = False, bottom=False)
plt.tight_layout()

#######################################################################

import torch
from model import *

# for small networks, CPU is faster
device = "cpu"

# current device for testing
print(f"Using {device} device")

# load and instrument trained model
trained = NeuralNetwork()
trained.load_state_dict(torch.load("model.pth"))
model = InstrumentedNetwork(trained, nodes).to(device)
model.eval()
print(model)

#######################################################################

from dataset import *

# test data and success bitmap
data = HandwrittenDigitsDataset('t10k')

# split test images according to label
digit = [data.image[data.label == i] for i in range(10)]

#######################################################################

import matplotlib.animation as animation

# called to advance animation to next frame
def animate(i):
	with torch.no_grad(): model(digit[8][i].unsqueeze(0).to(device))

animation = animation.FuncAnimation(fig, animate, frames=360, interval=1000/6)
#animation.save('fire-8.mp4')

# show or save the figure
plt.show() if file is None else plt.savefig(file, bbox_inches='tight', pad_inches=0.0, transparent=True)
