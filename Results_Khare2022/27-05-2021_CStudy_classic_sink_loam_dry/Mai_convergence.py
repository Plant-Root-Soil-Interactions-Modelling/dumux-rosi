
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from math import *
import csv
from pathlib import Path



#col = ["r--", "k:", "tab:olive", "tab:blue", "tab:pink", "tab:orange",  "tab:gray", "b", "m-", "c", "darkslategrey", "g", "k" ]

col = [ "k:","b" ,"tab:blue","g"]


col_WCT = ["w","w","r", "tab:olive", "tab:blue", "tab:pink", "tab:orange",  "tab:gray", "b", "m", "c", "darkslategrey", "g"]
#os.walk('cumT_plots/.')
#data = np.loadtxt("cumT_plots/PyBind_Loam_szenB_Tend3_dt120_10_4cm",delimiter=';')
#pl, = ax1.plot(data[0,:],data[1,:], color="r") # TPot



x = 3  # time for outputs [days]
scenario = "loam_dry"
x2= 2.99
y=  [0,        0,          3.4,          4.7,         5.625,          6.95,          8.5,            9.525,           10.55,        13.75,        15.3,      18.05,       19.4]
WCT=[str("0"),str("0"),str("1.9 min"),str("21.2 h"),str("6.4 h"),str("1.3 h"),str("6.8 min"),str("2.4 min"),str("1.4 min"),str("39 s"), str("28 s"),str("23 s"),str("19 s")]
# Txt-structure:   np.savetxt(dirname + "_pressures", np.vstack((x_, mean_bulk_soil_pressure, ESWP, cpx, max_bulk_soil_pressure, min_bulk_soil_pressure, mean_xylem_pressure)), delimiter = ';')
gridsize = ["4.0 cm", "2.0 cm", "1.0 cm", "reference"] #,  "1.0 cm Rhizo"]
gridsize_inverted = gridsize[::-1]
# Tact
l=[]
#discretization_error = ["", "", str("(0.30)"), str("(0.49)"), str("(0.62)"), str("(1.00)"), str("(1.24)"), str("(1.48)"), str("(2.24)"), str("(2.60)"), str("(3.26)"), "", ""]
fig, (ax1) = plt.subplots(figsize=(7,6))
#ax1.set_title('cumT Convergence Study Clay dry')
l=[]
#test = 2160*8*3-1 #constant timestep to get values at t=3
for dirname, dirnames, filenames in os.walk('Convergence_mai/.'): # this transpiration file is without skip, otherwise integration of transpiration over skip-time steps is wrong. skip needs to be there to get rid of the spike in the actual transpiration graph.
    filenames.sort()
    for i,f in enumerate(filenames):
        try:
            l.append(Path(f).stem)
            data = np.loadtxt("Convergence_mai/"+f,delimiter=';')  

            pl, = ax1.plot(data[0,:],data[3,:],col[i], label = gridsize_inverted[i]) #+" " +str(discretization_error[i]))
            #ax1.text(x2, y[i], WCT[i], fontsize=8, color=col_WCT[i], ha="right")
            #print(gridsize_inverted[i], (data[0,test],data[3,test])) # dt120
            #print(gridsize_inverted[i], (data[0,4318],data[3,4318])) # dt60
        except Exception as ex:
            print(ex)  
#pl, = ax1.plot(data[0,:],data[3,:],"k") # TPot draw (bad coding, but easy)

ax1.set_xlabel("Time [days]")  
ax1.set_ylabel("Cumulative transpiration [cm³]")  
ax1.legend( loc="best", ncol=2, fontsize=9)

ax1.set_xlim(0,x)
ax1.set_ylim(0, 4.5)
#ax1.set_title("(A)")
#to read discretization error
#ax1.set_yticks(minor_ticks)

#plt.tight_layout()
plt.savefig(scenario+"_CumT_Mai.jpg", dpi=300, bbox_inches = "tight")
plt.show()

