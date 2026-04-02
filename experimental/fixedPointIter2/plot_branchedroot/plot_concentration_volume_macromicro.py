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

font = {'size'   : 20}
plt.rc('font', **font)
mpl.rcParams['mathtext.default'] = 'regular'
path2file = "../scripts/results_branchedroot/Exudate_medium1/"

left  = 0.1  # the left side of the subplots of the figure
right = 0.85    # the right side of the subplots of the figure
bottom = 0.1   # the bottom of the subplots of the figure
top = 0.95      # the top of the subplots of the figure
wspace = 0.3   # the amount of width reserved for blank space between subplots
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

def get_axis_limits(ax, scale=0.9):
    return ax.get_xlim()[1]*scale, ax.get_ylim()[1]*scale


root_types = ['straight1', 'branched1']
titles = ['Straight root', 'Branched root']
annots = ['(a)','(b)']
figs = ['straight_root', 'branched_root']

fig, ax = plt.subplots(1,len(root_types))

for m in range(0, len(root_types)): 
        
    scenario = root_types[m]+'/'

    #macro 
    res = 1
    with open(path2file + scenario + "TotC_macro.csv") as f:
        totC = [list(map(float, row)) for row in csv.reader(f)]
    max_len = max(len(row) for row in totC)
    totC_ = np.array([row + [np.nan]*(max_len-len(row)) for row in totC])    
    conc_macro = totC_ / (res**3) 

    with open(path2file + scenario + "WaterC_macro.csv") as f:
        wc = [list(map(float, row)) for row in csv.reader(f)]
    wc_ = np.array([row + [np.nan]*(max_len-len(row)) for row in wc])  
    watvol_macro = wc_ * (res**3)
    time = np.loadtxt(path2file + scenario + "time.txt", delimiter=",")[:-1,0]


    #micro
    path_cyl = path2file+"/"+scenario+"cyl_val/"
    folder = Path(path_cyl)
    count = len(list(folder.glob("*time*.txt")))
    for i in range(0, count): 
        with open(path_cyl + 'Cyl_cellVol_'+str(i)+".txt") as f:
            cellvol_ = [list(map(float, row)) for row in csv.reader(f)]
        max_len = max(len(row) for row in cellvol_)
        cellvol = np.array([row + [np.nan]*(max_len-len(row)) for row in cellvol_])    
        
        with open(path_cyl + 'Cyl_content1_'+str(i)+".txt") as f:
            totC_ = [list(map(float, row)) for row in csv.reader(f)]
        totC = np.array([row + [np.nan]*(max_len-len(row)) for row in totC_])  
        conc_micro_ = np.divide(totC,cellvol)
        
        with open(path_cyl + 'Cyl_watercontent_'+str(i)+".txt") as f:
            wc_ = [list(map(float, row)) for row in csv.reader(f)]
        wc = np.array([row + [np.nan]*(max_len-len(row)) for row in wc_])  
        watvol_micro_ = np.multiply(wc,cellvol)
            
        if i == 0: 
            max_rows = np.shape(cellvol)[0]
            conc_micro = conc_micro_
            watvol_micro = watvol_micro_
        else: 
            conc_micro_ = pad_rows_top(conc_micro_, max_rows)
            conc_micro = np.hstack([conc_micro, conc_micro_])
            watvol_micro_ = pad_rows_top(watvol_micro_, max_rows)
            watvol_micro = np.hstack([watvol_micro, watvol_micro_])
            
    conc_ = np.hstack([conc_macro, conc_micro])
    watvol_ = np.hstack([watvol_macro, watvol_micro])
    idx = np.argsort(-conc_, axis = 1) 
    conc = np.take_along_axis(conc_, idx, axis = 1)
    watvol = np.take_along_axis(watvol_, idx, axis = 1)
    watvol_cum = np.cumsum(watvol, axis=1)

    colors = plt.cm.jet(np.linspace(0,1,np.shape(conc)[0]))
    for i in range(0, np.shape(conc)[0]): 
        if np.around(int(time[i] *1000)/1000-int(time[i]),2) == 0.5 : #only plot lunch time concentration/volumes for better visibility
            x = conc[i,:]
            y = watvol_cum[i,:]
            line = ax[m].plot(x,y, color = colors[i])
        # ax[1].set_title(add_labels[i][j])   
        

    # ax[m].set_xscale('log')
    # ax[m].set_yscale('log')
    ax[m].set_xlim([0,6e-5])
    ax[m].set_ylim([0, 15])
    # ax[m].set_title(titles[m])
    ax[m].annotate(annots[m], xy=get_axis_limits(ax[m]))
    
    #insert paraview screenshots
    arr_img = plt.imread('plots/'+figs[m]+".png")
    im = OffsetImage(arr_img, zoom=0.75, interpolation='auto')
    ab = AnnotationBbox(im, (1, 0), xycoords='axes fraction', box_alignment=(1.3,-0.03), frameon = False)
    ax[m].add_artist(ab)
            
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.5, 0.02, 0.4])
sm = plt.cm.ScalarMappable(cmap=cm.jet,norm=plt.Normalize(vmin=0, vmax=int(time[-1])))
clb = fig.colorbar(sm, cax=cbar_ax)
clb.ax.set_title('Time (d)\n')


fig.supxlabel('Exudate concentration ($mol/cm^3$ water)', y = 0.38)
fig.supylabel('Water volume ($cm^3$)', x = 0.07, y = 0.75)
# plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
plt.subplots_adjust(top=0.95, bottom=0.5)
plt.show()

