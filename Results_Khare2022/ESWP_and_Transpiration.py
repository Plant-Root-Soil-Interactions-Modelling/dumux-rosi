
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



x = 3   # time for outputs [days]
scenario = "loam_dry"

# Txt-structure:   np.savetxt(dirname + "_pressures", np.vstack((x_, mean_bulk_soil_pressure, ESWP, cpx, max_bulk_soil_pressure, min_bulk_soil_pressure, mean_xylem_pressure)), delimiter = ';')
gridsize = ["Potential transpiration" , "4.0 cm", "3.0 cm", "2.0 cm", "1.5 cm", "1.0 cm", "0.8 cm", "0.6 cm", "0.4 cm"] #,  "1.0 cm Rhizo"]
gridsize_inverted = gridsize[::-1]
# Tact
l=[]
for dirname, dirnames, filenames in os.walk('cumT_plots/.'):
    filenames.sort()
    for i,f in enumerate(filenames):
        try:
            l.append(Path(f).stem)
            data = np.loadtxt("cumT_plots/"+f,delimiter=';')
            #pl, = ax1.plot(data[0,:],data[1,:],col[i]) # TPot
            pl, = ax1.plot(data[0,:],data[1,:],col[i], label = gridsize_inverted[i]) # Tact
            
            
        except Exception as ex:
            print("Something went wrong with file "+f)    
            print(ex)
pl, = ax1.plot(data[0,:],data[2,:],"k") # TPot draw (bad coding, but easy)
ax1.set_xlabel("Time [days]")  
ax1.set_ylabel("Actual transpiration [cm³/day]")

plt.legend(bbox_to_anchor=(0,-0.16,1,0.2), loc="lower left",mode="expand", borderaxespad=0, ncol=4)

ax1.set_xlim(0,x)
#ax1.set_title('Tact Convergence Study Clay dry')

plt.tight_layout()
plt.savefig(scenario+"_Tact.png", dpi=300, bbox_inches = "tight")
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

            pl, = ax1.plot(data[0,:],data[3,:],col[i], label = gridsize_inverted[i])
            
        except Exception as ex:
            print(ex)  
#pl, = ax1.plot(data[0,:],data[3,:],"k") # TPot draw (bad coding, but easy)
ax1.set_xlabel("Time [days]")  
ax1.set_ylabel("Cumulative transpiration [cm³]")  
plt.legend(bbox_to_anchor=(0,-0.16,1,0.2), loc="lower left",mode="expand", borderaxespad=0, ncol=4)
ax1.set_xlim(0,x)
plt.tight_layout()
plt.savefig(scenario+"_CumT.png", dpi=300, bbox_inches = "tight")
plt.show()






# ESWP
l=[]
for dirname, dirnames, filenames in os.walk('ESWP/.'):
    filenames.sort()
    for i,f in enumerate(filenames):
        try:
            l.append(Path(f).stem)
            data = np.loadtxt("ESWP/"+f,delimiter=';')
            #pl, = ax1.plot(data[0,:],data[1,:],col[i]) # TPot
            #pl, = ax1.plot(data[0,:],data[2,:],col[i]) # ESWP
            plt.plot(data[0,:],data[2,:],col[i], label = gridsize_inverted[i], linewidth=1.5) # ESWP
            #pl, = ax1.plot(data[0,:],data[3,:],col[i]) # root collar pressure
            
        except Exception as ex:
            print("Something went wrong with file "+f)    
            print(ex)
plt.xlabel("Time [days]")  
plt.ylabel("ESWP [cm]")
gridsize = ["4.0 cm", "3.0 cm", "2.0 cm", "1.5 cm", "1.0 cm", "0.8 cm", "0.6 cm", "0.4 cm"] #,  "1.0 cm Rhizo"]
gridsize_inverted = gridsize[::-1]
plt.legend(bbox_to_anchor=(0,-0.28,1,0.2), loc="lower left",mode="expand", borderaxespad=0, ncol=4)
plt.xlim(0,x)
plt.ylim(-15000, 0)
#ax1.set_title('Tact Convergence Study Clay dry')
plt.tight_layout()
plt.savefig(scenario+"_ESWP_all.png", dpi=300, bbox_inches = "tight")
print(l)
plt.show()



# ESWP per grid

l=[]
for dirname, dirnames, filenames in os.walk('ESWP/.'):
    filenames.sort()
    for i,f in enumerate(filenames):
        try:
            l.append(Path(f).stem)
            data = np.loadtxt("ESWP/"+f,delimiter=';')
            #plt.plot(0,0, label = str(gridsize_inverted[i]) +":")
            plt.plot(data[0,:],data[1,:],col[i],linestyle ="solid", label = "mean bulk soil pressure") # mean bulk soil pressure
            plt.plot(data[0,:],data[2,:],col[i], linestyle ="dashed", label = "ESWP") # ESWP
            plt.plot(data[0,:],data[3,:],col[i],linestyle ="dotted", label = "root collar pressure") # root collar pressure
            plt.xlabel("Time [days]")  
            plt.ylabel("pressure head [cm]")
            gridsize = ["4.0 cm", "3.0 cm", "2.0 cm", "1.5 cm", "1.0 cm", "0.8 cm", "0.6 cm", "0.4 cm"] #,  "1.0 cm Rhizo"]
            gridsize_inverted = gridsize[::-1]
            #plt.legend(ncol = 5, prop ={"size":9}, loc='lower center')
            plt.legend(bbox_to_anchor=(0,-0.2,1,0.2), loc="lower left",mode="expand", borderaxespad=0, ncol=4)

            plt.xlim(0,x)
            plt.ylim(-15000, 0)
            #ax1.set_title('Tact Convergence Study Clay dry')
            plt.tight_layout()
            plt.savefig(scenario+"_ESWP_" + str(gridsize_inverted[i]) + "_same_color.png", dpi=300, bbox_inches = "tight")
            plt.show()
            


            plt.plot(data[0,:],data[1,:],"brown",linestyle ="solid", label = "bulk soil pressure") # mean bulk soil pressure
            plt.plot(data[0,:],data[2,:],"b", linestyle ="dashed", label = "ESWP") # ESWP
            plt.plot(data[0,:],data[3,:],"g",linestyle ="dotted", label = "root collar pressure") # root collar pressure
            sdfdfplt.plot((0,-15000),(3,-15000), linestyle = solid)
            plt.xlabel("Time [days]")  
            plt.ylabel("pressure head [cm]")
            gridsize = ["4.0 cm", "3.0 cm", "2.0 cm", "1.5 cm", "1.0 cm", "0.8 cm", "0.6 cm", "0.4 cm"] #,  "1.0 cm Rhizo"]
            gridsize_inverted = gridsize[::-1]
            plt.legend(bbox_to_anchor=(0,-0.2,1,0.2), loc="lower left",mode="expand", borderaxespad=0, ncol=4)

            plt.xlim(0,x)
            plt.ylim(-15000, 0)
            #ax1.set_title('Tact Convergence Study Clay dry')
            plt.tight_layout()
            plt.savefig(scenario+"_ESWP_" + str(gridsize_inverted[i]) + "_different_color.png", dpi=300, bbox_inches = "tight")
            plt.show()

        except Exception as ex:
            print("Something went wrong with file "+f)    
            print(ex)
