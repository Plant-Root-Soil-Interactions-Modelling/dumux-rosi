import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from scipy import interpolate
import pandas as pd
from os import walk
import re
import csv
from pathlib import Path
import os

fontsize = 24
mpl.rcParams.update({
    'font.size': fontsize,
    'axes.titlesize': fontsize,
    'axes.labelsize': fontsize,
    'xtick.labelsize': fontsize,
    'ytick.labelsize': fontsize,
    'legend.fontsize': fontsize,
    'figure.titlesize': fontsize,
    'mathtext.default': 'regular'
})

comp = ['flowonly', 'numeric']
cols = [['b', 'r'],['c', 'm']]
linw = [1,3]
fig, ax = plt.subplots(1)

for i in range(0, len(comp)): 

    path2file = "../scripts/results_"+comp[i]+"/"

    soiltype = ['loam', 'sand']
    diffusion = ['low', 'medium', 'mediumhigh', 'high']
    sorption = ['low', 'medium', 'high']
    

    
    scenario = soiltype[0]+'/'

    #time
    time = np.array(np.loadtxt(path2file + scenario + "time.txt", delimiter=",")[:-1,0])
    #theta
    theta = np.array(np.loadtxt(path2file + scenario + 'theta.csv', delimiter=',').mean(axis=1))
    mintheta = np.array(np.loadtxt(path2file + scenario + 'theta.csv', delimiter=',').min(axis=1))
    maxtheta = np.array(np.loadtxt(path2file + scenario + 'theta.csv', delimiter=',').max(axis=1))
    #pHead
    pHead = np.array(np.loadtxt(path2file + scenario + 'pHead.csv', delimiter=',').mean(axis=1))
    minpHead = np.array(np.loadtxt(path2file + scenario + 'pHead.csv', delimiter=',').min(axis=1))
    maxpHead = np.array(np.loadtxt(path2file + scenario + 'pHead.csv', delimiter=',').max(axis=1))


    ax.plot(time, theta[:len(time)],color = cols[i][0],  linewidth = linw[i])
    ax.plot(time, mintheta[:len(time)],color = cols[i][0],linestyle = '--',  linewidth = linw[i])
    ax.plot(time, maxtheta[:len(time)],color = cols[i][0],linestyle = ':',  linewidth = linw[i])
    # ax.set_xlim(left = 22)
    ax2 = ax.twinx()
    ax2.plot(time, pHead[:len(time)], color = cols[i][1], linewidth = linw[i])
    ax2.plot(time, minpHead[:len(time)],color = cols[i][1],linestyle = '--',  linewidth = linw[i
    ])
    ax2.plot(time, maxpHead[:len(time)],color = cols[i][1],linestyle = ':',  linewidth = linw[i])
    ax2.set_ylim([-4000,0])

# fig.suptitle(titles[0])
ax.plot(np.nan, np.nan, 'k-', label = 'Mean') 
ax.plot(np.nan, np.nan, 'k--', label = 'Min') 
ax.plot(np.nan, np.nan, 'k:', label = 'Max') 
ax.plot(np.nan, np.nan, color = cols[0][0], label = 'theta, analytic') 
ax.plot(np.nan, np.nan, color = cols[0][1], label = 'pHead, analytic') 
ax.plot(np.nan, np.nan, color = cols[1][0], label = 'theta, numeric') 
ax.plot(np.nan, np.nan, color = cols[1][1], label = 'pHead, numeric') 

ax.set_xlabel('Time (d)', x = 0.47, y = 0.03, fontsize=fontsize)
ax.set_ylabel('Soil water content ($cm^3$ $cm^{-3}$)',  fontsize=fontsize)
ax2.set_ylabel('Soil  water potential ($cm$)',  fontsize=fontsize)
# ax2.legend(loc = 0,bbox_to_anchor=(1,0.9))
ax.legend(loc = 'lower left')


plt.show()

