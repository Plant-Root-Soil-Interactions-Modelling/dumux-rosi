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

path = '../../scripts/results/Exudate/loam_2_exu_nodecay_withsorption/swp_arrays/'
files = [f for f in listdir(path) if isfile(join(path, f))]
files = np.sort(files) 
print(files) 

fig, axs = plt.subplots(1, len(files))
for i in range(0, len(files)): 
    conc = np.load(path + files[i])
    # conc_ = conc[:,int(np.shape(conc)[1]/2),600:]
    conc_ = np.max(conc, axis = 1) 
    conc_ = conc_[:,600:]
    conc_rot = ndimage.rotate(conc_, 90)
    colors = plt.cm.jet(conc_)
    if i == 0: 
        norm = mpl.colors.Normalize(vmin=np.min(conc_), vmax=np.max(conc_))
    m = axs[i].imshow(conc_rot, cmap = 'jet', norm = norm, aspect=1)
    if i == len(files)-1: 
        divider = make_axes_locatable(axs[i])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(m, cax = cax, label = 'Exudate concentration (-)')
    # plt.title("Exudate concentration")
    axs[i].axis('off')
plt.show()