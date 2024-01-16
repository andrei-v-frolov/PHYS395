#!/usr/bin/env python
# introduction to plotting in Python
#######################################################################

from os import environ
import numpy as np
from pltconfig import *
import matplotlib.pyplot as plt

# show figure in macOS window
#environ['matplotlib.backend'] = 'macosx'

# output figure as PDF file
environ['matplotlib.backend'] = 'pdf'

#######################################################################

n = 32

x = np.linspace(0.0,10.0,n)
y = np.sin(x)

print(x)
print(y)

plt.fill_between(x,y,'c', alpha=0.1)
plt.plot(x,y,'r.-')

plt.xlim(x[0],x[-1])

#plt.show()
plt.savefig("output.pdf", bbox_inches='tight', pad_inches=0.02, transparent=True)
