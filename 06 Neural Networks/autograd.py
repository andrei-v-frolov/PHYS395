#!/usr/bin/env python
# PyTorch compute graph and back-propagation demo

#######################################################################

import torch

#######################################################################

x = torch.ones(5)  # input tensor
y = torch.zeros(3)  # expected output
w = torch.randn(5, 3, requires_grad=True)
b = torch.randn(3, requires_grad=True)
z = torch.matmul(x, w)+b
loss = torch.nn.functional.binary_cross_entropy_with_logits(z, y)

print(f"Gradient function for z = {z.grad_fn}")
print(f"Gradient function for loss = {loss.grad_fn}")

# backpropagation
loss.backward()
print(w.grad)
print(b.grad)

# disabling gradient
print(z.requires_grad)

with torch.no_grad():
    z1 = torch.matmul(x, w)+b
print(z1.requires_grad)

z2 = z.detach()
print(z2.requires_grad)