# matplotlib configuration for Planck papers
# adopted from https://github.com/zonca/paperplots

import os
import numpy as np
from matplotlib import rcParams, rc

# unit conversion
def cm2inch(cm):
    """Centimeters to inches"""
    return cm/2.54

# do you really want TeX?
def usetex(flag=True):
    """Shortcut to enable/disbleable TeX typesetting"""
    rcParams.update({'text.usetex': flag})

# common setup for matplotlib
params = {'backend': os.environ.get('matplotlib.backend', 'pdf'),
          'savefig.dpi': 300, # save figures to 300 dpi
          'axes.labelsize': 16,
          'legend.fontsize': 16,
          'xtick.labelsize': 16,
          'ytick.major.pad': 11,
          'xtick.major.pad': 11,
          'ytick.labelsize': 16,
          'text.usetex': True,
          'font.family':'sans-serif',
          'font.size': 16,
          # free font similar to Helvetica
          'font.sans-serif':'FreeSans'}

# use of Sans Serif also in math mode
rc('text.latex', preamble=r'\usepackage{sfmath}')

rcParams.update(params)

# import pyplot by default
import matplotlib.pyplot as plt
