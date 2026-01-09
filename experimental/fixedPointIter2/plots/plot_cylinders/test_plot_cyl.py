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

path = '../../scripts/results/Exudate/loam_4_exu_nodecay_nosorption/exu_arrays/'
conc = np.load(path + 'day55.npy')

# conc_ = conc[:,int(np.shape(conc)[1]/2),:]
conc_ = np.max(conc, axis = 1) 
conc_rot = ndimage.rotate(conc_, 90)
colors = plt.cm.seismic_r(conc_)
norm = mpl.colors.Normalize(vmin=np.min(conc_), vmax=np.max(conc_))
plt.imshow(conc_rot, cmap = 'seismic_r', norm = norm, aspect=1)#1.67
m = cm.ScalarMappable(cmap=plt.cm.seismic_r, norm=norm)
#m.set_array([])
#plt.colorbar(m, ax = ax2, label = 'Water content (-)')
plt.title("Exudate concentration")
plt.axis('off')
plt.show()