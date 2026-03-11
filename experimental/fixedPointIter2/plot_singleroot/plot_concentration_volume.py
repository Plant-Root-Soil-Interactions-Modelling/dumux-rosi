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

font = {'size'   : 16}
plt.rc('font', **font)
mpl.rcParams['mathtext.default'] = 'regular'
path2file2 = "../scripts/results_singleroot/Exudate"

left  = 0.1  # the left side of the subplots of the figure
right = 0.85    # the right side of the subplots of the figure
bottom = 0.1   # the bottom of the subplots of the figure
top = 0.95      # the top of the subplots of the figure
wspace = 0.3   # the amount of width reserved for blank space between subplots
hspace = 0.5   # the amount of height reserved for white space between subplots

def get_num(filename): 
    regex = re.compile(r'\d+')
    return [int(x) for x in regex.findall(filename)]

#Scenarios 
soil1 =['loam', 'loam','loam', 'loam']
soil2 =['sand', 'loam','loam', 'loam']
sorption1 = ['low', 'low', 'low', 'low']
sorption2 = ['low', 'medium', 'low', 'low']
swp1 =['100','100','100', '100'] 
swp2 =['100','100', '100', '1000'] 
trans1 = ['0', '0', '0', '0.6']
trans2 = ['0', '0', '0.6', '0.6']
titles = ['Soil type', 'Sorption', 'Transpiration', 'Initial soil water potential']
add_labels = [['Loam', 'Sand'],['Low sorption', 'Medium sorption'], ['Transpiration = 0 cm/d', 'Transpiration = 0.6 cm/d'], ['Initial SWP = -100 cm', 'Initial SWP = -1000 cm']]
dt = 20/60/24 #day

fig, ax = plt.subplots(len(soil1),2)
for i in range(0, len(soil1)): 
    scenario1 = soil1[i]+'_res1_sorption'+sorption1[i]+'_SWP_ini'+swp1[i]+'_trans'+trans1[i]
    scenario2 = soil2[i]+'_res1_sorption'+sorption2[i]+'_SWP_ini'+swp2[i]+'_trans'+trans2[i]
    
    path1 = path2file2+"/"+scenario1+"/exu_arrays/"
    path2 = path2file2+"/"+scenario2+"/exu_arrays/"
    filenames = next(walk(path1), (None, None, []))[2]
    filenames = sorted(filenames, key = get_num)
    colors = plt.cm.jet(np.linspace(0,1,len(filenames)))
    
    for m in range(0, len(filenames)): 
        df1 = np.load(path1+filenames[m])
        df2 = np.load(path2+filenames[m])
        df = [df1, df2]
        
        #exudate balance input - output
        for j in range(0,2): 
            x = -np.sort(-df[j].flatten())
            y = np.arange(1, len(x)+1) / 1000
            line = ax[i,j].plot(x,y, color = colors[m])
            ax[i,j].set_title(add_labels[i][j])   
            ax[i,j].set_xscale('log')
            ax[i,j].set_yscale('log')
            # ax[i,j].invert_xaxis()
            ax[i,j].set_xlim([1e-8, 1e-3])
            # ax[i,j].set_ylim([0,300])
            
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.9, 0.15, 0.02, 0.7])
sm = plt.cm.ScalarMappable(cmap=cm.jet,norm=plt.Normalize(vmin=0, vmax=int(float(get_num(filenames[-1])[0])/10)))
clb = fig.colorbar(sm, cax=cbar_ax)
clb.ax.set_title('Time (d)')



fig.supxlabel('Exudate concentration ($mol/cm^3$)', y = 0.01)
fig.supylabel('Soil volume ($cm^3$)', x = 0.03)
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
plt.show()

