#!/usr/bin/env python
# plot stock prices for trading days
#######################################################################

import numpy as np
import matplotlib.pyplot as plt

stocks = np.loadtxt('stocks.dat')
n,columns = stocks.shape

plt.plot(stocks)
plt.legend(["AAPL", "AMD", "AMZN", "CSCO", "GOOGL", "INTC", "META", "MSFT", "NFLX", "NVDA", "QCOM", "SBUX", "TSLA"])
plt.xlim([0,n])

plt.show()
