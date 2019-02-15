# matplotlib configuration for Planck papers
# adopted from https://github.com/zonca/paperplots

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
params = {'backend': 'pdf',
          'savefig.dpi': 300, # save figures to 300 dpi
          'axes.labelsize': 10,
          'legend.fontsize': 10,
          'xtick.labelsize': 10,
          'ytick.major.pad': 6,
          'xtick.major.pad': 6,
          'ytick.labelsize': 10,
          'text.usetex': True,
          'font.family':'sans-serif',
          'font.size': 10,
          # free font similar to Helvetica
          'font.sans-serif':'FreeSans'}

# use of Sans Serif also in math mode
rc('text.latex', preamble='\usepackage{sfmath}')

rcParams.update(params)

# import pyplot by default
import matplotlib.pyplot as plt
