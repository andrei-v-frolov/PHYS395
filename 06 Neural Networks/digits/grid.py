# addressable computed grid of pixel coordinates

#######################################################################

import numpy as np

#######################################################################

class grid:
	def __init__(self, shape, dx=(1,0), dy=(0,1), center=(0,0), **kwargs):
		ny,nx = shape
		self.nx = nx
		self.ny = ny
		self.index = 0
		self.count = nx*ny
		self.shape = (ny,nx)
		self.dx = np.array(dx)
		self.dy = np.array(dy)
		self.center = np.array(center)
		self.kwargs = kwargs
	
	# indexing runs in [y,x] order (as in plt.imshow)
	def __getitem__(self, idx):
		i,j = idx
		x =  (j-(self.nx-1)/2)*self.dx
		y = -(i-(self.ny-1)/2)*self.dy
		return self.center + y + x
	
	# restart iterator counter
	def __iter__(self):
		self.index = 0
		return self
	
	# iterator varies x fastes (as in np.flatten)
	def __next__(self):
		if self.index >= self.count:
			raise StopIteration
		i = self.index//self.nx
		j = self.index - i*self.nx
		self.index += 1
		return self[i,j]
	
	# pixel coordinates as NumPy array
	def array(self):
		return np.fromiter(self, dtype=np.dtype((float,2)))
