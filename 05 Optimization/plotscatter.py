#!/usr/bin/env python
# Plot MCMC-sampled probability density function
# usage: plotscatter [width list (in cm)]

###############################################################################
# import libraries
###############################################################################

# configure import path
import os, sys, warnings

# import libraries
from math import *
import numpy as np
from scipy import stats
from scipy.interpolate import interp1d
from sys import argv, stdin

# plot style options
aspect = 4/4.   # figure aspect ratio (default is 4:3)
fill = False    # produce transparent plots if false
grid = False    # do we want to render the plot grid?

widths = map(float, sys.argv[1].split(',')) if (len(argv) > 1) else [18.0, 12.0, 8.8]

###############################################################################
# import sampled data
###############################################################################

data = np.loadtxt('DATA')

x = data[:,0]
y = data[:,1]

###############################################################################
# create the plots
###############################################################################

# Configure Matplotlib options for pretty output
from pltconfig import *
from matplotlib import gridspec
from matplotlib.ticker import MaxNLocator

def plot_setup(ax, xlabel='', ylabel='', xlim=[2,20], ylim=[-4,0]):
    """Setup common plot style parameters"""
    
    # labels
    if (xlabel):
        plt.xlabel(xlabel)
    if (ylabel):
        plt.ylabel(ylabel)
    
    # distance of axis label to tick labels
    ax.yaxis.labelpad = 10*width/17.
    ax.xaxis.labelpad = 10*width/17.
    
    # reduce ticks for small figures
    if width < 10:
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
    # grid
    if (grid):
        plt.grid(True, which="major", axis="both")
    
    # axes limits
    plt.ylim(ylim)
    plt.xlim(xlim)
    
    # set vertical y axis ticklables
    for ticklabel in ax.yaxis.get_ticklabels():
        ticklabel.set_rotation("vertical")

# Create the plots
for width in widths:
    fig = plt.figure(figsize=(cm2inch(width), cm2inch(width/aspect)), frameon=fill)
    
    ax = fig.add_subplot(111)
    
    #ax.plot(x, y, 'k.', markersize=1, zorder=1)
    ax.plot(x[::3], y[::3], 'k.', markersize=1, zorder=1)
    
    xmin,xmax = plt.xlim(); ymin,ymax = plt.ylim()
    xmin,xmax = [-2,2]; ymin,ymax = [-2,2]
    X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    
    # kernel density estimator
    positions = np.vstack([X.ravel(), Y.ravel()])
    scatter = np.vstack([x, y])
    kernel = stats.gaussian_kde(scatter)
    Z = np.reshape(kernel(positions).T, X.shape)
    
    # confidence level lookup
    Q = np.sort(Z.flatten()); S = np.cumsum(Q)
    P = interp1d(1.0-S/S[-1], Q, kind='linear')
    
    ax.imshow(np.rot90(Z), cmap=plt.cm.Blues, alpha=0.8, extent=[xmin, xmax, ymin, ymax], aspect='auto', zorder=0)
    ax.contour(X, Y, Z, [P(0.95),P(0.68)], colors=('b'), zorder=5)
    
    #plt.title("distribution of Gaussian peak CDF parameters")
    plot_setup(ax, xlabel=r'$x$', ylabel=r'$y$', xlim=[xmin,xmax], ylim=[ymin,ymax])
    
    # reduce white space around figure
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0, hspace=0)
    
    # save to pdf with right bounding box
    base = "mcmc-pdf"
    pdf = "%s.%dmm.pdf" % (base,int(width*10)) if len(widths) > 1 else base+".pdf"
    plt.savefig(pdf, bbox_inches='tight', pad_inches=0.02, transparent=not(fill))
