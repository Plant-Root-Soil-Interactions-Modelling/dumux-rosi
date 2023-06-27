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
path0 = "../results/concentration/"
path1 = "../results/"
data0 = np.load(path0+"/day028.npy", allow_pickle=True)
data1 = np.load(path1+"/soilc_maize_exudate_2019.npy", allow_pickle=True)
#############################################################
params = {'legend.fontsize': 18,
          'xtick.labelsize': 18,
          'ytick.labelsize': 18,
          'legend.handlelength': 2}

plt.rcParams.update(params)

widtha = 20
widthb = 44
depth = 74

res = 0.1

nx1 = int(widtha / res);
ny1 = int(widthb / res);
nz1 = int(depth / res);

nx2 = 10;
ny2 = 22;
nz2 = 37;


cyl = np.reshape(data0,(nx1,ny1,nz1))
soil = np.reshape(data1[-1,:],(nz2,ny2,nx2))

fig, (ax1, ax2) = plt.subplots(1,2)

A = (np.amax(cyl,axis=0));
A = np.rot90(A)
A = np.flipud(A)
im1 = ax1.imshow(A,  cmap=cm.jet,
               origin='lower', extent=[0, widtha, 0, depth],
                 vmax=abs(A).max(), vmin=10**-30)
divider2 = make_axes_locatable(ax1)
cax1 = divider2.append_axes("right", size="10%", pad=.1)
cbar1 = plt.colorbar(im1, cax=cax1, format="%.2E")
cbar1.set_label('Exudate concentration at day 6 (g/cm3)', rotation=270,size=18, labelpad=30)
ax1.xaxis.set_visible(False)
ax1.yaxis.set_visible(False)
ax1.annotate('(a)', xy=get_axis_limits(ax1),color = 'w',size=18)


C = (np.amax(soil,axis=2));
#C = np.rot90(C)
#C = np.flipud(C)
C = C.astype(float)
im2 = ax2.imshow(C,  cmap=cm.jet,
               origin='lower', extent=[0, widtha, 0, depth],
                 vmax=abs(C).max(), vmin=10**-30)
divider1 = make_axes_locatable(ax2)
cax2 = divider1.append_axes("right", size="10%", pad=0.05)
cbar2 = plt.colorbar(im2, cax=cax2, format="%.2E")
cbar2.set_label('Exudate concentration at day 8 (g/cm3)', rotation=270,size=18, labelpad=40)
ax2.xaxis.set_visible(False)
ax2.yaxis.set_visible(False)
ax2.annotate('(b)', xy=get_axis_limits(ax2),color = 'w',size=18)


plt.show()

