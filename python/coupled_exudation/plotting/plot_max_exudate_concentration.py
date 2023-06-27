import sys; 
sys.path.append("../../../..")
import os
import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import matplotlib.colors as colors
import math
#import plantbox as rb
import re
import os
import scipy.sparse
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm
from matplotlib.ticker import MultipleLocator
from matplotlib.pylab import *
#############################################################
def get_axis_limits(ax, scale=.1):
    return ax.get_xlim()[1]*0.8, ax.get_ylim()[1]*0.1

def get_numbers_from_filename(filename):
    return re.search(r'\d+', filename).group(0)
#############################################################
#load rhizodeposit volumes 
path = "../results/concentration/"
data0 = np.load(path+"/day006.npy")
data1 = np.load(path+"/day008.npy")
data2 = np.load(path+"/day010.npy")
#############################################################
params = {'legend.fontsize': 18,
          'xtick.labelsize': 18,
          'ytick.labelsize': 18,
          'legend.handlelength': 2}

plt.rcParams.update(params)

widtha =20
widthb = 44
depth = 74

res = 0.1
nx = int(widtha / res);
ny = int(widthb / res);
nz = int(depth / res);
cols = ['b', 'g']

day006 = np.reshape(data0,(nx,ny,nz))
day008 = np.reshape(data0,(nx,ny,nz))
day010 = np.reshape(data0,(nx,ny,nz))

fig, (ax1, ax2, ax3) = plt.subplots(1,3)

A = (np.amax(day006,axis=0));
A = np.rot90(A)
A = np.flipud(A)
im1 = ax1.imshow(A, interpolation='bilinear', cmap=cm.jet,
               origin='lower', extent=[0, widtha, 0, depth],
               vmax=abs(A).max(), vmin=10**-15)
divider2 = make_axes_locatable(ax1)
cax1 = divider2.append_axes("right", size="10%", pad=.1)
cbar1 = plt.colorbar(im1, cax=cax1, format="%.0f")
cbar1.set_label('Exudate concentration at day 6 (g/cm3)', rotation=270,size=18, labelpad=30)
ax1.xaxis.set_visible(False)
ax1.yaxis.set_visible(False)
ax1.annotate('(a)', xy=get_axis_limits(ax1),color = 'w',size=18)



C = (np.amax(day008,axis=0));
C = np.rot90(C)
C = np.flipud(C)
im2 = ax2.imshow(C, interpolation='bilinear', cmap=cm.jet,
               origin='lower', extent=[0, widtha, 0, depth],
                 vmax=abs(C).max(), vmin=10**-15)
divider1 = make_axes_locatable(ax2)
cax2 = divider1.append_axes("right", size="10%", pad=0.05)
cbar2 = plt.colorbar(im2, cax=cax2, format="%.0f")
cbar2.set_label('Exudate concentration at day 8 (g/cm3)', rotation=270,size=18, labelpad=40)
ax2.xaxis.set_visible(False)
ax2.yaxis.set_visible(False)
ax2.annotate('(b)', xy=get_axis_limits(ax2),color = 'w',size=18)

M = (np.amax(day010,axis=0));
M = np.rot90(M)
M = np.flipud(M)
im3 = ax3.imshow(M, interpolation='bilinear', cmap=cm.jet,
               origin='lower', extent=[0, widtha, 0, depth],
               vmax=abs(M).max(), vmin=10**-15)
divider2 = make_axes_locatable(ax3)
cax3 = divider2.append_axes("right", size="10%", pad=.1)
cbar3 = plt.colorbar(im3, cax=cax3, format="%.0f")
cbar3.set_label('Exudate concentration at day 10 (g/cm3)', rotation=270,size=18, labelpad=30)
ax3.xaxis.set_visible(False)
ax3.yaxis.set_visible(False)
ax3.annotate('(c)', xy=get_axis_limits(ax3),color = 'w',size=18)

plt.show()

