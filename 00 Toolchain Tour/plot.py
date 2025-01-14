#!/usr/bin/env python
# introduction to plotting in Python

#######################################################################

import numpy as np

# interval to sample on
n = 32; a = 0.0; b = 10.0

x = np.linspace(a, b, n)
y = np.sin(x)

print(x)
print(y)

#######################################################################

import matplotlib.pyplot as plt

plt.figure()

plt.plot(x, y, 'o-')
plt.axhline(0.0, color='black', linestyle=':', zorder=0)
plt.fill_between(x, y, -1.0, color='green', alpha=0.2)

plt.xlim([a,b])

plt.show()
