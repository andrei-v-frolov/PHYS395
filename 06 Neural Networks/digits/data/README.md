# [MNIST database of handwritten digits](http://yann.lecun.com/exdb/mnist/)

### Supplied data

- `train-images-idx3-ubyte.gz` - 60,000 training images, IDX-3, gzipped
- `train-labels-idx1-ubyte.gz` - 60,000 training labels, IDX-1, gzipped
- `t10k-images-idx3-ubyte.gz` - 10,000 test images, IDX-3, gzipped
- `t10k-labels-idx1-ubyte.gz` - 10,000 test labels, IDX-1, gzipped

### IDX-1 data format

~~~
[offset] [type]          [value]          [description]
0000     32 bit integer  0x00000801(2049) magic number
0004     32 bit integer  60000            number of labels
0008     unsigned byte   ??               label
0009     unsigned byte   ??               label
........
xxxx     unsigned byte   ??               label
~~~

### IDX-3 data format

~~~
[offset] [type]          [value]          [description]
0000     32 bit integer  0x00000803(2051) magic number
0004     32 bit integer  60000            number of images
0008     32 bit integer  28               number of rows
0012     32 bit integer  28               number of columns
0016     unsigned byte   ??               pixel
0017     unsigned byte   ??               pixel
........
xxxx     unsigned byte   ??               pixel
~~~
