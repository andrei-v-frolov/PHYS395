# Recognition of handwritten digits

### Data products

- `data` - [MNIST database of handwritten digits](http://yann.lecun.com/exdb/mnist/)
- `models` - saved models (after a few epochs of training)
- `images` - training, test, and success rate renders

### Neural network implementation

- `idx.py` - idx data format parser
- `plot.py` - render training and test MNIST handwritten data images
- `dataset.py` - PyTorch dataset provider for MNIST handwritten data
- `model.py` - PyTorch neural network model for MNIST handwritten data
- `train.py` - train, save, and test handwritten digit recognition network
- `test.py` - render success rate on test MNIST handwritten data images
