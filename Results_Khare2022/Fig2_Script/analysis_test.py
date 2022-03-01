
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from math import *
import csv
from pathlib import Path
import matplotlib.gridspec as gridspec


# Figure 2 (7 days)
fig = plt.figure(figsize=(13,6))
spec = gridspec.GridSpec(ncols=3, nrows=1, figure=fig)

ax1 = fig.add_subplot(spec[0, 0])
ax2 = fig.add_subplot(spec[0, 1])
ax3 = fig.add_subplot(spec[0, 2])
col = ["r--", "tab:pink", "tab:orange",  "tab:gray", "b", "m-", "c", "darkslategrey", "g", "g-"]

# 
l=[]
for dirname, dirnames, filenames in os.walk('cumT_plots_sand_dry/.'):
    filenames.sort()
    for i,f in enumerate(filenames):
        try:
            l.append(Path(f).stem)
            data = np.loadtxt('cumT_plots_sand_dry/'+f,delimiter=';')
            ax1.plot(data[0,:],data[3,:],col[i])
            
            
        except Exception as ex:
            print("Something went wrong with file "+f)    
            print(ex)
            
# loam
l=[]
for dirname, dirnames, filenames in os.walk('cumT_plots_loam_dry/.'):
    filenames.sort()
    for i,f in enumerate(filenames):
        try:
            l.append(Path(f).stem)
            data = np.loadtxt('cumT_plots_loam_dry/'+f,delimiter=';')
            ax2.plot(data[0,:],data[3,:],col[i])
            
            
        except Exception as ex:
            print("Something went wrong with file "+f)    
            print(ex)
            
# clay
l=[]
for dirname, dirnames, filenames in os.walk('cumT_plots_clay_dry/.'):
    filenames.sort()
    for i,f in enumerate(filenames):
        try:
            l.append(Path(f).stem)
            data = np.loadtxt('cumT_plots_clay_dry/'+f,delimiter=';')
            ax3.plot(data[0,:],data[3,:],col[i])
            
            
        except Exception as ex:
            print("Something went wrong with file "+f)    
            print(ex)
            
ax1.set_xlabel("Time [days]")  
ax1.set_ylabel("Cumulative transpiration [cm³]")
ax2.set_xlabel("Time [days]")  
ax2.set_ylabel("Cumulative transpiration [cm³]")
ax3.set_xlabel("Time [days]")  
ax3.set_ylabel("Cumulative transpiration [cm³]")

gridsize = ["4.0 cm", "3.0 cm", "2.0 cm", "1.5 cm", "1.0 cm", "0.8 cm", "0.6 cm", "0.4 cm",  "1.0 cm multiscale"]
gridsize_inverted = gridsize[::-1]

fig.legend(gridsize_inverted, bbox_to_anchor=(0.51,-0.05), loc="lower center", ncol=9, prop={'size':12}) 
#plt.subplots_adjust(bottom=0.02)
plt.tight_layout()

ax1.set_title("(A)")
ax2.set_title("(B)")
ax3.set_title("(C)")
plt.tight_layout()

ax1.set_xlim(0,7)
ax1.set_ylim(0,4)

ax2.set_xlim(0,7)
ax2.set_ylim(0,50)

ax3.set_xlim(0,7)
ax3.set_ylim(0,50)

plt.savefig('Figure_2.jpg', dpi=300, bbox_inches = "tight")
plt.show()
