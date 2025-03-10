# PyTorch model for MNIST handwritten data

#######################################################################

# import PyTorch libraries
import torch
from torch import nn

#######################################################################

# get CPU, GPU or MPS device for training
device = (
    "cuda" if torch.cuda.is_available()
    else "mps" if torch.backends.mps.is_available()
    else "cpu"
)

# neural network implementation
class NeuralNetwork(nn.Module):
    def __init__(self):
        super().__init__()
        self.flatten = nn.Flatten()
        self.linear_sigmoid_stack = nn.Sequential(
            nn.Linear(28*28, 11*11),
            nn.Sigmoid(),
            nn.Linear(11*11, 5*5),
            nn.Sigmoid(),
            nn.Linear(5*5, 10)
        )

    def forward(self, x):
        x = self.flatten(x)
        logits = self.linear_sigmoid_stack(x)
        return logits
