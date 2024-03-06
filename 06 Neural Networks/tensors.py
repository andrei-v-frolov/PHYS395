#!/usr/bin/env python
# PyTorch tensor operations demo

#######################################################################

import torch
import numpy as np

#######################################################################

# initialize from data
data = [[1, 2],[3, 4]]
x = torch.tensor(data)
print(x)

# initialize from numpy array
array = np.array(data)
x = torch.from_numpy(array)
print(x)

# tensor initializers
ones = torch.ones_like(x)
print(f"Ones Tensor: \n {ones} \n")

rand = torch.rand_like(x, dtype=torch.float)
print(f"Random Tensor: \n {rand} \n")

# explicit shape initializers
shape = (2,3,)
rand = torch.rand(shape)
ones = torch.ones(shape)
zero = torch.zeros(shape)

print(f"Random Tensor: \n {rand} \n")
print(f"Ones Tensor: \n {ones} \n")
print(f"Zeros Tensor: \n {zero}")

# tensor attributes
tensor = torch.rand(3,4)

print(f"Shape of tensor: {tensor.shape}")
print(f"Datatype of tensor: {tensor.dtype}")
print(f"Device tensor is stored on: {tensor.device}")

# get CPU, GPU or MPS device
device = (
    "cuda" if torch.cuda.is_available()
    else "mps" if torch.backends.mps.is_available()
    else "cpu"
)
print(f"\nUsing {device} device")

# move tensor to GPU (if available)
tensor = tensor.to(device)

print(f"Shape of tensor: {tensor.shape}")
print(f"Datatype of tensor: {tensor.dtype}")
print(f"Device tensor is stored on: {tensor.device}")

# indexing and slicing
tensor = torch.ones(3, 4) # .to(device)
print(f"First row: {tensor[0]}")
print(f"First column: {tensor[:, 0]}")
print(f"Last column: {tensor[..., -1]}")
tensor[:,1] = 0
print(tensor)

# concatenation
t1 = torch.cat([tensor, tensor, tensor], dim=1)
print(t1)

# flatten tensor
print(tensor)
f = tensor.flatten()
print(f)

# matrix multiplication between two tensors (y1, y2, y3 will have the same value)
# tensor.T returns the transpose of a tensor
y1 = tensor @ tensor.T
y2 = tensor.matmul(tensor.T)
y3 = torch.rand_like(y1)
torch.matmul(tensor, tensor.T, out=y3)
print(y1)

# element-wise product (z1, z2, z3 will have the same value)
z1 = tensor * tensor
z2 = tensor.mul(tensor)
z3 = torch.rand_like(tensor)
torch.mul(tensor, tensor, out=z3)
print(z1)

# single-element tensor
agg = tensor.sum()
item = agg.item()
print(item, type(item))

# in-place operations
print(f"{tensor} \n")
tensor.add_(5)
print(tensor)

# bridge to numpy
t = torch.ones(5)
print(f"t: {t}")
n = t.numpy()
print(f"n: {n}")

# bridged array traces changes
t.add_(1)
print(f"t: {t}")
print(f"n: {n}")
