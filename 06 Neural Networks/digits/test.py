#!/usr/bin/env python
# render success rate on test MNIST handwritten data images

# parse arguments
from sys import argv
file = argv[1] if len(argv) > 1 else None

# import PyTorch libraries
import torch
from model import *
from dataset import *

# current device for testing
print(f"Using {device} device")

# test data
test_data = HandwrittenDigitsDataset('t10k')
test_dataloader = DataLoader(test_data, batch_size=1, shuffle=False, pin_memory=True)

# load trained model
model = NeuralNetwork().to(device)
model.load_state_dict(torch.load("model.pth"))
model.eval()
print(model)

# success bitmap
data = idx.open('data/t10k-images-idx3-ubyte.gz')
success = torch.zeros((len(test_data),28,28))

with torch.no_grad():
	for i, (X, y) in enumerate(test_dataloader):
		X, y = X.to(device), y.to(device)
		pred = model(X)
		success[i,:,:] = 1.0 if (pred.argmax(1) == y).item() else 0.0

# import matplotlib
import matplotlib
matplotlib.use('macosx' if file is None else 'PDF')

from pylab import *

image = idx.grid(data, 125, 80)
success = idx.grid(success, 125, 80)

# create the figure
figure(figsize=(10,8), frameon=False)
s = matplotlib.colors.LinearSegmentedColormap.from_list("success", ["red", "white"])
d = matplotlib.colors.LinearSegmentedColormap.from_list("success", ["#00000000", "black"])
imshow(success, origin='upper', cmap=s, vmin=0.0, vmax=1.0, interpolation='none')
imshow(image, origin='upper', cmap=d, vmin=0.0, vmax=255.0, interpolation='none')

tick_params(left=False, right=False, labelleft=False, labelbottom = False, bottom=False)

if not(file is None): plt.savefig(file, bbox_inches='tight', pad_inches=0.0, transparent=True)
show()