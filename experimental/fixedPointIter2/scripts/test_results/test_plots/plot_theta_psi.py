import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib
matplotlib.use("QtAgg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from scipy import interpolate
import pandas as pd
from os import walk
import re
import csv
from pathlib import Path

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
path2file = "../"

soiltype = ['loam', 'sand']
diffusion = ['low', 'medium', 'mediumhigh', 'high']
sorption = ['low', 'medium', 'high']
cols = ['b', 'r', 'g']

fig, ax = plt.subplots(1)
scenario = soiltype[1]+'_diffusion'+diffusion[2]+'_sorption'+sorption[0]+'/'

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
#psiXyl
psiXyl = []
minpsiXyl = []
maxpsiXyl = []
with open(path2file + scenario+'psiXyl.txt', "r") as f:
    for line in f:
        values = [float(x) for x in line.split(",")]
        if len(values) == 0:
            continue
        psiXyl.append(np.mean(values))
        minpsiXyl.append(np.min(values))
        maxpsiXyl.append(np.max(values))
psiXyl = np.array(psiXyl)
minpsiXyl = np.array(minpsiXyl)
maxpsiXyl = np.array(maxpsiXyl)

ax.plot(time, theta[:len(time)],color = cols[0],  label = 'Mean soil water content')
ax.plot(time, mintheta[:len(time)],color = cols[0],linestyle = '--',  label = 'Min soil water content')
ax.plot(time, maxtheta[:len(time)],color = cols[0],linestyle = ':',  label = 'Max soil water content')
ax2 = ax.twinx()
ax2.plot(time, pHead[:len(time)], color = cols[1],  label = 'Mean pHead')
ax2.plot(time, minpHead[:len(time)],color = cols[1],linestyle = '--',  label = 'Min pHead')
ax2.plot(time, maxpHead[:len(time)],color = cols[1],linestyle = ':',  label = 'Max pHead')

ax2.plot(time, psiXyl[:len(time)], color = cols[2],  label = 'Mean psiXyl')
ax2.plot(time, minpsiXyl[:len(time)],color = cols[2],linestyle = '--',  label = 'Min psiXyl')
ax2.plot(time, maxpsiXyl[:len(time)],color = cols[2],linestyle = ':',  label = 'Max psiXyl')
ax2.set_ylim([-1000,10])


# fig.suptitle(titles[0])
ax.set_xlabel('Time (d)', x = 0.47, y = 0.03, fontsize=fontsize)
ax.set_ylabel('Soil water content ($cm^3$ $cm^{-3}$)',  fontsize=fontsize)
ax2.set_ylabel('Soil / xylem water potential ($cm$)',  fontsize=fontsize)
# ax2.legend(loc = 0,bbox_to_anchor=(1,0.9))
ax.legend(loc = 'lower left')
ax2.legend(loc = 'upper left')

plt.show()

