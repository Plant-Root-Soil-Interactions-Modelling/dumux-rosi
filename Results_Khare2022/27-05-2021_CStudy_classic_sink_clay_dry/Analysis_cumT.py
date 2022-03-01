
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from math import *
import csv
from pathlib import Path

fig, ax1 = plt.subplots(figsize=(8,7)) 
col = ["tab:pink", "tab:orange",  "tab:gray", "b", "m-", "c", "y", "g", "g-"]
#os.walk('cumT_plots/.')
#data = np.loadtxt("cumT_plots/PyBind_Loam_szenB_Tend3_dt120_10_4cm",delimiter=';')
#pl, = ax1.plot(data[0,:],data[1,:], color="r") # TPot

# Tact
l=[]
for dirname, dirnames, filenames in os.walk('cumT_plots/.'):
    filenames.sort()
    for i,f in enumerate(filenames):
        try:
            l.append(Path(f).stem)
            data = np.loadtxt("cumT_plots/"+f,delimiter=';')
            #pl, = ax1.plot(data[0,:],data[1,:],col[i]) # TPot
            pl, = ax1.plot(data[0,:],data[1,:],col[i]) # Tact
            
            
        except Exception as ex:
            print("Something went wrong with file "+f)    
            print(ex)
pl, = ax1.plot(data[0,:],data[2,:],"k") # TPot draw (bad coding, but easy)
ax1.set_xlabel("Time [days]")  
ax1.set_ylabel("Actual transpiration [cm³/day]")
gridsize = ["Potential transpiration" , "4.0 cm", "3.0 cm", "2.0 cm", "1.5 cm", "1.0 cm", "0.8 cm", "0.6 cm", "0.4 cm"] #,  "1.0 cm Rhizo"]
gridsize_inverted = gridsize[::-1]
ax1.legend(gridsize_inverted, loc="upper right") 

ax1.set_xlim(0,14)
#ax1.set_title('Tact Convergence Study Clay dry')

plt.tight_layout()
plt.savefig('clay_dry_Tact.png', dpi=300, bbox_inches = "tight")
print(l)
#plt.show()

fig, ax1 = plt.subplots(figsize=(6,7)) 
#ax1.set_title('cumT Convergence Study Clay dry')
l=[]

for dirname, dirnames, filenames in os.walk('cumT_plots/.'):
    filenames.sort()
    for i,f in enumerate(filenames):
        try:
            l.append(Path(f).stem)
            data = np.loadtxt("cumT_plots/"+f,delimiter=';')  

            pl, = ax1.plot(data[0,:],data[3,:],col[i])
            
        except Exception as ex:
            print(ex)  
#pl, = ax1.plot(data[0,:],data[3,:],"k") # TPot draw (bad coding, but easy)
ax1.set_xlabel("Time [days]")  
ax1.set_ylabel("Cumulative transpiration [cm³]")  
ax1.legend(gridsize_inverted, loc="upper left") 
ax1.set_xlim(0,14)


plt.tight_layout()
plt.savefig('clay_dry_CumT.png', dpi=300, bbox_inches = "tight")
plt.show()
