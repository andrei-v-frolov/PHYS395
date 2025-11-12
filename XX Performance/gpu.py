#!/usr/bin/env python
# time finite difference Laplacian stencils on GPU using PyTorch

#######################################################################

import torch
import numpy as np

# get CPU, GPU or MPS device
device = (
    "cuda" if torch.cuda.is_available()
    else "mps" if torch.backends.mps.is_available()
    else "cpu"
)
print(f"\nUsing '{device}' device...")

from time import perf_counter as now

#######################################################################

fields = 1; n = 1000000

# test data for stencil application
x = torch.zeros([fields,n], device=device); x[:,3] = 1.0

# second order Laplacian stencil in 1D
stencil = torch.nn.Conv1d(fields, fields, kernel_size=3, bias=False, device=device)
stencil.weight.data[0,0] = torch.tensor([1,-2,1])

with torch.no_grad():
    t1 = now(); y = stencil(x); t2 = now()
    print(f"1D PyTorch convolve on {y.device}", y[0,0:7].cpu().numpy(), t2-t1, sep='\n')

#######################################################################

n = 1000

# test data for stencil application
x = torch.zeros([fields,n,n], device=device); x[:,3,3] = 1.0

# second order Laplacian stencil in 2D
stencil = torch.nn.Conv2d(fields, fields, kernel_size=3, bias=False, device=device)
stencil.weight.data[0,0] = torch.tensor([[0,1,0],[1,-4,1],[0,1,0]])

# second order isotropic Laplacian stencil in 2D
#stencil = np.array([[1,4,1],[4,-20,4],[1,4,1]])/6.0

with torch.no_grad():
    t1 = now(); y = stencil(x); t2 = now()
    print(f"2D PyTorch convolve on {y.device}", y[0,0:7,0:7].cpu().numpy(), t2-t1, sep='\n')


#######################################################################
