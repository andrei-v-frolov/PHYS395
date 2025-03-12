#!/usr/bin/env python
# train, save, and test handwritten digit recognition network

#######################################################################

# import PyTorch libraries
import torch
from model import *
from dataset import *

#######################################################################

# for small networks, CPU is faster
device = "cpu"

# current device for training
print(f"Using {device} device")

#######################################################################

# stochastic batch size
batch = 25

# training and test data
train = HandwrittenDigitsDataset('train')
test = HandwrittenDigitsDataset('t10k')

# data loader for training
loader = DataLoader(train, batch_size=batch, shuffle=True, pin_memory=True)

#######################################################################

# model definition
model = NeuralNetwork().to(device)
print(model)

# model parameters optimizer
loss_fn = nn.CrossEntropyLoss()
optimizer = torch.optim.SGD(model.parameters(), lr=0.75)
scheduler = torch.optim.lr_scheduler.ExponentialLR(optimizer, gamma=0.9)

# train network for a single epoch
def train(loader, model, loss_fn, optimizer, log_every=100):
    count = len(loader.dataset); model.train()
    for batch, (X, y) in enumerate(loader):
        X, y = X.to(device), y.to(device)
        
        # compute prediction error
        pred = model(X)
        loss = loss_fn(pred, y)
        
        # backpropagation
        loss.backward()
        optimizer.step()
        optimizer.zero_grad()
        
        # log progress every few batches
        if (batch + 1) % log_every == 0:
            loss, current = loss.item(), (batch + 1) * len(X)
            print(f"loss: {loss:>7f}  [{current:>5d}/{count:>5d}]")

# check the model performance against the test dataset
def validate(data, model, loss_fn):
    count = len(data); model.eval()
    with torch.no_grad():
        pred = model(data.image.to(device))
        loss = loss_fn(pred, data.label.to(device))
        success = (pred.argmax(-1).cpu() == data.label).sum()
    print(f"Success rate: {(100*success/count):>0.1f}%, avg loss: {loss:>8f} \n")

#######################################################################

# training process
epochs = 10
for t in range(epochs):
    print(f"Epoch {t+1}\n-------------------------------")
    train(loader, model, loss_fn, optimizer)
    validate(test, model, loss_fn)
    torch.save(model.state_dict(), f"model-{t}.pth")
    scheduler.step()
print("Done!")

# save trained model
torch.save(model.state_dict(), "model.pth")
print("Saved PyTorch model state to model.pth")
