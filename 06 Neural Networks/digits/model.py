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
        self.stack = nn.Sequential(
            nn.Flatten(),
            nn.Linear(28*28, 11*11),
            nn.Sigmoid(),
            nn.Linear(11*11, 5*5),
            nn.Sigmoid(),
            nn.Linear(5*5, 10)
        )

    def forward(self, x):
        logits = self.stack(x)
        return logits
