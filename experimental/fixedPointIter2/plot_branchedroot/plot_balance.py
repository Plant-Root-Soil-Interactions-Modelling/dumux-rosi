import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from scipy import interpolate
import pandas as pd
from os import walk
import re
import csv
from pathlib import Path
from matplotlib.offsetbox import OffsetImage, AnnotationBbox

font = {'size'   : 16}
plt.rc('font', **font)
path2file = "../scripts/results_branchedroot/Exudate/"


titles = ['Straight root', 'Branched root']
annots = [['(a)','(b)'],['(c)','(d)'],['(e)','(f)']]
figs = ['straight_root', 'branched_root']

linst = ['-', '--', ':', '-.']
cols= ['r', 'b', 'g', 'y']
linst = ['-', '--']
dt = 20/60/24 #day

fig, ax = plt.subplots(1)

scenario1 = 'straight/'
scenario2 = 'branched/'
df1 = pd.read_csv(path2file+"/"+scenario1+"/exud.csv")
df2 = pd.read_csv(path2file+"/"+scenario2+"/exud.csv")
df = [df1, df2]

#exudate balance input - output
x_sim = np.arange(0, dt*len(df1['Exud_tot'].loc[:].values),dt)  
ax.plot(x_sim, df1['Exud_tot'].loc[:].values,color = cols[0])
ax.plot(x_sim, df1['Exud_ads'].loc[:].values,  color = cols[1])
ax.plot(x_sim,  df1['Exud_liq'].loc[:].values, color = cols[2])
ax.plot(x_sim, df1['Exud_decay'].loc[:].values, color = cols[3])

ax.plot(x_sim, df2['Exud_tot'].loc[:].values,linestyle = '--',color = cols[0])
ax.plot(x_sim, df2['Exud_ads'].loc[:].values,linestyle = '--',color = cols[1])
ax.plot(x_sim, df2['Exud_liq'].loc[:].values, linestyle = '--',color = cols[2])
ax.plot(x_sim, df2['Exud_decay'].loc[:].values,linestyle = '--',color = cols[3])


lines = ax.get_lines()
legend1 = ax.legend([lines[i] for i in np.arange(0,4)], ['Total amount of exudates', 'Sorbed exudates','Dissolved exudates', 'Decomposed exudates'], loc='upper right')
ax.add_artist(legend1)

for j in range(0, 2): 
    ax.plot([], [], color = 'k', linestyle = linst[j], label = titles[j])
ax.legend(loc = 'upper left')



fig.supxlabel('Time (d)')
fig.supylabel('Exudate content in soil domain (mol)', x = 0.05)
plt.show()

