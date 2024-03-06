# PyTorch dataset provider for MNIST handwritten data

# import PyTorch libraries
import torch
from torch.utils.data import Dataset
from torch.utils.data import DataLoader

# import IDX reader library
import idx
import numpy as np

# MNIST handwritten digits dataset
class HandwrittenDigitsDataset(Dataset):
	def __init__(self, set):
		self.image = torch.from_numpy(idx.open('data/' + set + '-images-idx3-ubyte.gz')/np.single(255.0))
		self.label = torch.from_numpy(idx.open('data/' + set + '-labels-idx1-ubyte.gz'))
		assert self.image.shape[0] == self.label.shape[0], "Dataset size mismatch..."

	def __len__(self):
		return self.image.shape[0]

	def __getitem__(self, idx):
		return self.image[idx], self.label[idx]
