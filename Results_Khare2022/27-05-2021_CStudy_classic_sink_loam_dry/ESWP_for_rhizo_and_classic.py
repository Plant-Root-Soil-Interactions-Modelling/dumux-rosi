
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from math import *
import csv
from pathlib import Path

fig, ax1 = plt.subplots(figsize=(8,7)) 
col = ["b", "r"]
#os.walk('cumT_plots/.')
#data = np.loadtxt("cumT_plots/PyBind_Loam_szenB_Tend3_dt120_10_4cm",delimiter=';')
#pl, = ax1.plot(data[0,:],data[1,:], color="r") # TPot



x = 3   # time for outputs [days]
scenario = "loam_dry"

# Txt-structure:   np.savetxt(dirname + "_pressures", np.vstack((x_, mean_bulk_soil_pressure, ESWP, cpx, max_bulk_soil_pressure, min_bulk_soil_pressure, mean_xylem_pressure)), delimiter = ';')
gridsize = ["multiscale", "1.0 cm"] #,  "1.0 cm Rhizo"]
gridsize_inverted = gridsize[::-1]


# ESWP per grid

l=[]
for dirname, dirnames, filenames in os.walk('ESWP_rhizo_classic/.'):
    filenames.sort()
    for i,f in enumerate(filenames):
        try:
            l.append(Path(f).stem)
            data = np.loadtxt("ESWP_rhizo_classic/"+f,delimiter=';')
            #plt.plot(0,0, label = str(gridsize_inverted[i]) +":")
            plt.plot(data[0,:],data[1,:],col[i],linestyle ="solid", label = "mean bulk soil pressure") # mean bulk soil pressure
            plt.plot(data[0,:],data[2,:],col[i], linestyle ="dashed",label = "ESWP") # ESWP
            #plt.plot(data[0,:],data[3,:],col[i],linestyle ="solid", label = "root collar pressure") # root collar pressure
            plt.xlabel("Time [days]")  
            plt.ylabel("pressure head [cm]")
            gridsize = ["1.0 cm", "multiscale"] #,  "1.0 cm Rhizo"]
            gridsize_inverted = gridsize[::-1]
            #plt.legend(ncol = 5, prop ={"size":9}, loc='lower center')

        except Exception as ex:
            print("Something went wrong with file "+f)    
            print(ex)
    plt.legend(bbox_to_anchor=(0,-0.135,1,0.135), loc="lower left",mode="expand", borderaxespad=0, ncol=4)

    plt.xlim(0,x)
    plt.ylim(-15000, 0)
    #ax1.set_title('Tact Convergence Study Clay dry')plt.tight_layout()
    plt.savefig(scenario+"_ESWP__rhizo_classic" + str(gridsize_inverted[i]) + "_same_color.png", dpi=300, bbox_inches = "tight")
    plt.show()
            
