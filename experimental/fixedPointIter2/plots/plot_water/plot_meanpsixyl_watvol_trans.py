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
import ast

font = {'size'   : 18}
plt.rc('font', **font)

fig, axs = plt.subplots(3,1)

# path = '../../scripts/results/Exudate/day1_loam_2_exu_withdecay_withsorption_strangesorption/'
path = '../../scripts/results/Exudate/day1_loam_2_noexu_day10/'
time = np.array(([x.split(',')[0] for x in open(path + 'time.txt').readlines()])).astype(float)[1:]

#Transpiration 
trans = np.array(([x.split(' ')[0] for x in open(path + 'trans.txt').readlines()])).astype(float)
axs[0].plot(time, trans, 'Actual transpiration') 
#potential transpiration 
axs[0].legend()
axs[0].set_title('Cumulative transpiration (cm^3)')


#Xylem potential 
mean_psiXyl = []
with open(path + "psiXyl.txt", "r") as f:
    for line in f:
        mean_psiXyl.append(np.mean(np.fromstring(line, sep=",")))
axs[1].plot(time, mean_psiXyl) 
axs[1].set_title('Mean xylem water potential (cm)')

#Water volume 
watVol_raw = np.loadtxt(path + "watVol.csv", delimiter=",")
domainvol = np.shape(watVol_raw)[1]*2**3
watVol = np.sum(watVol_raw, axis = 1)
thetas = 0.494
axs[2].plot(time, watVol, label = 'water volume in soil domain') 
axs[2].plot(time, thetas*domainvol*np.ones((len(time))), label = 'Full saturation') 
axs[2].legend()
axs[2].set_title('Total water volume in soil domain (cm^3)')


plt.show()
