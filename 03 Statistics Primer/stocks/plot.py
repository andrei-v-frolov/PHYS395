#!/usr/bin/env python
# introduction to plotting in Python
#######################################################################

import numpy as np
import matplotlib.pyplot as plt

stocks = np.loadtxt('stocks.dat')

plt.plot(stocks)
plt.xlim([0,2516])

plt.show()
