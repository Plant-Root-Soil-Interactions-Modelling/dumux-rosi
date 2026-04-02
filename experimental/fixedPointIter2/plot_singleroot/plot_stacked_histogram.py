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

fontsize = 20
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
path2file = "../scripts/results_singleroot/Exudate/"

left  = 0.1  # the left side of the subplots of the figure
right = 0.85    # the right side of the subplots of the figure
bottom = 0.1   # the bottom of the subplots of the figure
top = 0.95      # the top of the subplots of the figure
wspace = 0.15   # the amount of width reserved for blank space between subplots
hspace = 0.5   # the amount of height reserved for white space between subplots

def get_num(filename): 
    regex = re.compile(r'\d+')
    return [int(x) for x in regex.findall(filename)]
    
def pad_rows_top(arr, target_rows):
    if arr.shape[0] < target_rows:
        # Create a NaN array of missing rows
        padding = np.full((target_rows - arr.shape[0], arr.shape[1]), np.nan)
        return np.vstack([padding, arr])  # NaNs on top
    return arr

def get_axis_limits(ax, scalex=0.1, scaley=0.8):
    return ax.get_xlim()[1]*scalex, ax.get_ylim()[1]*scaley


soiltype = ['loam', 'sand']
sorption = ['low', 'medium', 'high']
titles = ['Loam', 'Sand']
titles_y = ['Carboxylate sorption', 'Sugar sorption', 'Amino acid sorption']
annots = [['(a)','(b)'],['(c)','(d)'],['(e)','(f)']]
dt = 20/60/24 #day

fig, ax = plt.subplots(len(sorption),len(soiltype))

for m in range(0, len(soiltype)): 
    for n in range(0, len(sorption)): 
        
        scenario = soiltype[m]+'_res1_sorption'+sorption[n]+'_SWP_ini100_trans0/'
        df = pd.read_csv(path2file + scenario+"exud.csv")
        total = df['Exud_tot'].loc[:].values
        total[total==0] = 1e-10

        a = df['Exud_liq'].loc[:].values/total*100
        b = df['Exud_decay'].loc[:].values/total*100
        c = df['Exud_ads'].loc[:].values/total*100
        a[np.isnan(a)] = 100
        b[np.isnan(b)] = 0
        c[np.isnan(c)] = 0
        
        # print(a[-1], 'dissolved', soiltype[m], sorption[n])
        # print(b[-1], 'decomposed', soiltype[m], sorption[n])
        print(c[-1], 'sorbed', soiltype[m], sorption[n])
        
        data = np.vstack((a,b,c)).T
        index_ = np.arange(0, dt*len(df['Exud_tot'].loc[:].values),dt) 
        t_max = np.ceil(index_[-1])
        df_fin = pd.DataFrame(data,index=index_, columns=['dissolved exudates', 'decomposed exudates', 'sorbed exudates'])
        
        ax[n,m].set_xticks([0, t_max/2, t_max])
        ax[n,m].set_xlim(0, t_max)
        ax[n,m].set_ylim(0, 100)

        df_fin.plot.area(ax=ax[n,m])   
        ax[n,m].axvline(x=9, color='red', linestyle='--', linewidth=1.5)
        
        if n == 0 and m == 0:
            ax[n,m].legend(loc='upper right')
        else:
            ax[n,m].get_legend().remove()
        if n == 0: 
            ax[n,m].set_title(titles[m])
        if m == 0: 
            ax[n,m].set_ylabel(titles_y[n])
        ax[n,m].annotate(annots[n][m], xy=get_axis_limits(ax[n,m]))


fig.supxlabel('Time (d)', x = 0.47, y = 0.02, fontsize=fontsize)
fig.supylabel('Share of total root exudates (%)',  x = 0.015, fontsize=fontsize)
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
plt.show()

