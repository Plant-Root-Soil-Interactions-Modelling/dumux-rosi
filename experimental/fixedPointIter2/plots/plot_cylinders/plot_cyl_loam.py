import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import pandas as pd
from scipy import ndimage
import sys
from os import listdir
from os.path import isfile, join
from mpl_toolkits.axes_grid1 import make_axes_locatable

font = {'size'   : 18}
plt.rc('font', **font)

path = '../../scripts/results/Exudate/day1_loam_2_exu_withdecay_withsorption/'
arrays = ['swp_arrays/', 'exu_arrays/', 'decay_arrays/']
labels = ['Soil water potential (cm)', 'Exudate concentration (mol/cm^3)','Cumulative decay at a given location (mol/cm^3)']
files = [f for f in listdir(path+arrays[0]) if isfile(join(path+arrays[0], f))]
files = np.sort(files) 
 
fig, axs = plt.subplots(len(arrays), len(files))
for i in range(len(arrays)): 
    for j in range(0, len(files)): 
        data = np.load(path + arrays[i]+files[j])
        data_ = data[600:,int(np.shape(data)[1]/2),:]
        # data_ = np.max(data, axis = 1) 
        # data_ = data_[:,600:]
        data_rot = ndimage.rotate(data_, 180)
        colors = plt.cm.jet(data_)
        if j == 0 : 
            norm = mpl.colors.Normalize(vmin=np.min(data_), vmax=np.max(data_))
        # elif j ==0 & i == 1: 
            # norm = mpl.colors.Normalize(vmin=np.min(data_), vmax=np.max(data_))
        m = axs[i,j].imshow(data_rot, cmap = 'jet', norm = norm, aspect=1)
        if j == len(files)-1: 
            divider = make_axes_locatable(axs[i,j])
            cax = divider.append_axes("right", size="5%", pad=0.05)
            plt.colorbar(m, cax = cax, label = labels[i])
        axs[i,j].axis('off')
plt.show()