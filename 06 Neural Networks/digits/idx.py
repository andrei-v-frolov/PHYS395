# idx data format (for MNIST handwritten data)
# http://yann.lecun.com/exdb/mnist/

# load libraries
import gzip
import numpy as np

# load 32-bit unsigned integer from bytes (in MSB-first order)
def uint32_msb(data):
	return (data[0]<<24) + (data[1]<<16) + (data[2]<<8) + data[3]

# open gzipped file containing IDX-formatted uint8 array
def open(file):
	with gzip.open(file) as f:
		data = f.read()
		
		if not((data[0] == 0) and (data[1] == 0) and (data[2] == 8)):
			raise Exception("Invalid IDX magic number")
		
		dims = data[3]; shape = []
		for i in range(dims):
			shape.append(uint32_msb(data[4+(i<<2):8+(i<<2)]))
		
		return np.fromiter(data[4*(dims+1):], np.uint8).reshape(shape)

# make an image grid of data samples
def grid(array, nx, ny):
	(n,y,x) = array.shape
	assert nx*ny == n, "Grid shape does not conform to array size..."
	data = array.reshape([ny,nx,y,x]); image = np.zeros([ny*y,nx*x])

	for i in range(ny):
		for j in range(nx):
			image[i*y:(i+1)*y,j*x:(j+1)*x] = data[i][j]
	
	return image
