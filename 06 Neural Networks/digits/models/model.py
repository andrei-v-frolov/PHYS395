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
            nn.Linear(28*28, 11*11, bias=False),
            nn.Sigmoid(),
            nn.Linear(11*11, 5*5, bias=False),
            nn.Sigmoid(),
            nn.Linear(5*5, 10, bias=False)
        )

    def forward(self, x):
        return self.stack(x)

#######################################################################

# path-through module updating plot data
class Plot(nn.Module):
    def __init__(self, plot) -> None:
        super().__init__()
        self.plot = plot

    def forward(self, x):
        self.plot.set_array(x.cpu().numpy().flatten())
        return x

# network instrumented for plotting
class InstrumentedNetwork(NeuralNetwork):
    def __init__(self, trained, layer):
        super().__init__()
        self.stack = nn.Sequential(
            trained.stack[0],
            Plot(layer[0]),
            trained.stack[1],
            trained.stack[2],
            Plot(layer[1]),
            trained.stack[3],
            trained.stack[4],
            Plot(layer[2]),
            trained.stack[5],
            nn.Softmax(),
            Plot(layer[3])
        )

    def forward(self, x):
        return self.stack(x)
