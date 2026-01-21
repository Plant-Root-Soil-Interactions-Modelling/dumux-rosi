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

path = '../../scripts/results/Exudate/day1_loam_2_noexu/swp_arrays/'
files = [f for f in listdir(path) if isfile(join(path, f))]
files = np.sort(files) 
print(files) 

fig, axs = plt.subplots(1, len(files))
for i in range(0, len(files)): 
    swp = np.load(path + files[i])
    swp_ = swp[:,int(np.shape(swp)[1]/2),:]
    # swp_ = np.max(swp, axis = 1) 
    # swp_ = swp_[:,600:]
    swp_rot = ndimage.rotate(swp_, 180)
    colors = plt.cm.jet(swp_)
    if i == 0: 
        norm = mpl.colors.Normalize(vmin=np.min(swp_), vmax=np.max(swp_))
    m = axs[i].imshow(swp_rot, cmap = 'jet', norm = norm, aspect=1)
    if i == len(files)-1: 
        divider = make_axes_locatable(axs[i])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(m, cax = cax, label = 'Soil water potential (cm)')
    # plt.title("Exudate swpentration")
    axs[i].axis('off')
plt.show()