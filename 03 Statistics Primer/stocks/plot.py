#!/usr/bin/env python
# introduction to plotting in Python
#######################################################################

import numpy as np
import matplotlib.pyplot as plt

stocks = np.loadtxt('stocks.dat')

plt.plot(stocks)
plt.legend(["AAPL", "AMD", "AMZN", "CSCO", "GOOGL", "INTC", "META", "MSFT", "NFLX", "NVDA", "QCOM", "SBUX", "TSLA"])
plt.xlim([0,2516])

plt.show()
