
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from math import *
import csv
from pathlib import Path


import matplotlib.pyplot as plt
from numpy import genfromtxt
import numpy as np

data_a_act = np.genfromtxt("a_actual.out")
data_p = np.genfromtxt("a_pw_min_root.out")

collar_p = (data_p[:,1] - 1.e5) * 100. / 1.e3 / 9.81
ESWP_reference =  data_a_act[:,1]/0.0027 + collar_p

import os
import subprocess
import multiprocessing as mp
import matplotlib.pyplot as plt
from numpy import genfromtxt
import numpy as np

#data_tpot = np.genfromtxt("a_pot.out")
#data_a_act = np.genfromtxt("a_actual.out")
data_a_cumul = np.genfromtxt("a_cumul.out")
#data_b_act = np.genfromtxt("b_actual.out")
#data_b_cumul = np.genfromtxt("b_cumul.out")

#plot_color = ['r', 'g', 'b', 'r', 'g', 'b']

t = data_a_cumul[:,0]
dt = np.insert((t[1:] - t[:-1]), 0, 0)
#print(dt)
Alpha = 0.04
N = 1.6
Qr = 0.08
Qs = 0.43
M = 1-1/N

initialP = -652.3
NR = (1 + (Alpha * abs(initialP))**N)**M
WC_i = Qr + (Qs - Qr)/NR

print("Initial water content: ", WC_i)

DomainWater = WC_i * 960
print("Amount of water in domain: ", DomainWater)

WC = (DomainWater - data_a_cumul[:,1])/960
#print(WC)

head = 1/Alpha*(((Qs-Qr)/(WC - Qr))**(1/M) - 1)**(1/N)

fig, ax = plt.subplots()
ax.plot(dt, -head)
#plt.show()


# prepare plot
fig, ax = plt.subplots()
ax.plot(data_a_act[:,0], ESWP_reference, "k--", linewidth=1, label="reference $\psi_{s, eq}$")
ax.set_xlabel('Time [days]')
ax.set_ylabel('Water potential [cm]')
ax.legend(loc = 'upper right')
fig.savefig("eswp_reference.pdf", dpi=300, bbox_inches='tight')
plt.show()


col = ["r--", "k:", "tab:olive", "tab:blue", "tab:pink", "tab:orange",  "tab:gray", "b", "m-", "c", "darkslategrey", "g", "k" ]
col_WCT = ["w","w","r", "tab:olive", "tab:blue", "tab:pink", "tab:orange",  "tab:gray", "b", "m", "c", "darkslategrey", "g"]
#os.walk('cumT_plots/.')
#data = np.loadtxt("cumT_plots/PyBind_Loam_szenB_Tend3_dt120_10_4cm",delimiter=';')
#pl, = ax1.plot(data[0,:],data[1,:], color="r") # TPot



x = 3  # time for outputs [days]
scenario = "loam_dry"
x2= 2.99
y=  [0,        0,          3.4,          4.7,         5.625,          6.95,          8.5,            9.525,           10.55,        13.75,        15.3,      18.05,       19.4]
WCT=[str("0"),str("0"),str("4.3 min"),str("21.2 h"),str("6.4 h"),str("1.3 h"),str("6.8 min"),str("2.4 min"),str("1.4 min"),str("39 s"), str("28 s"),str("23 s"),str("19 s")]
# Txt-structure:   np.savetxt(dirname + "_pressures", np.vstack((x_, mean_bulk_soil_pressure, ESWP, cpx, max_bulk_soil_pressure, min_bulk_soil_pressure, mean_xylem_pressure)), delimiter = ';')
gridsize = [  "T$_{pot}$", "4.0 cm", "3.0 cm", "2.0 cm", "1.5 cm", "1.0 cm", "0.8 cm", "0.6 cm", "0.4 cm", "0.3 cm","0.2 cm", "reference" ,"1.0 cm multiscale"] #,  "1.0 cm Rhizo"]
gridsize_inverted = gridsize[::-1]
# Tact
l=[]

discretization_error = [str("(-0.20)"), "", str("(0.30)"), str("(0.49)"), str("(0.62)"), str("(1.00)"), str("(1.24)"), str("(1.48)"), str("(2.24)"), str("(2.60)"), str("(3.26)"), str("(3.57)"),""]
#test= 1
#test2= [0.0679235,  0.08118471, 0.09703501, 0.11597989, 0.13862351, 0.16568802, 0.19803653, 0.23670069, 0.28291354, 0.33814887]            
#plt.plot(test2, test)  

fig, (ax1,ax2) = plt.subplots(1,2, figsize=(13,6)) 
#ax1.set_title('cumT Convergence Study Clay dry')
l=[]
#test = 2160*8*3-1 #constant timestep to get values at t=3
for dirname, dirnames, filenames in os.walk('cumT_plots2/.'): # this transpiration file is without skip, otherwise integration of transpiration over skip-time steps is wrong. skip needs to be there to get rid of the spike in the actual transpiration graph.
    filenames.sort()
    for i,f in enumerate(filenames):
        try:
            l.append(Path(f).stem)
            data = np.loadtxt("cumT_plots2/"+f,delimiter=';')  

            pl, = ax1.plot(data[0,:],data[3,:],col[i], label = gridsize_inverted[i] +" " +str(discretization_error[i]))
            ax1.text(x2, y[i], WCT[i], fontsize=8, color=col_WCT[i], ha="right")
            #print(gridsize_inverted[i], (data[0,test],data[3,test])) # dt120
            #print(gridsize_inverted[i], (data[0,4318],data[3,4318])) # dt60
        except Exception as ex:
            print(ex)  
#pl, = ax1.plot(data[0,:],data[3,:],"k") # TPot draw (bad coding, but easy)

ax1.set_xlabel("Time [days]")  
ax1.set_ylabel("Cumulative transpiration [cm³]")  
ax1.legend( loc="best", ncol=2, fontsize=8)

ax1.set_xlim(0,x)
ax1.set_ylim(0, 20)
ax1.set_title("(A)")
#to read discretization error
#ax1.set_yticks(minor_ticks)

#plt.tight_layout()
plt.savefig(scenario+"_CumT.jpg", dpi=300, bbox_inches = "tight")


#fig, ax2 = plt.subplots(2,2) 

#os.walk('cumT_plots/.')
#data = np.loadtxt("cumT_plots/PyBind_Loam_szenB_Tend3_dt120_10_4cm",delimiter=';')
#pl, = ax1.plot(data[0,:],data[1,:], color="r") # TPot


col = ["b", "r"]
x = 3   # time for outputs [days]
scenario = "loam_dry"

# Txt-structure:   np.savetxt(dirname + "_pressures", np.vstack((x_, mean_bulk_soil_pressure, ESWP, cpx, max_bulk_soil_pressure, min_bulk_soil_pressure, mean_xylem_pressure)), delimiter = ';')
gridsize = ["multiscale", "1.0 cm"] #,  "1.0 cm Rhizo"]
gridsize_inverted = gridsize[::-1]


# ESWP per grid

l2=[]
for dirname, dirnames, filenames in os.walk('ESWP_rhizo_classic/.'):
    filenames.sort()
    for i,f in enumerate(filenames):
        try:
            l2.append(Path(f).stem)
            names = [ "1.0 cm RSS ","1.0 cm multiscale"]
            data = np.loadtxt("ESWP_rhizo_classic/"+f,delimiter=';')
            #plt.plot(0,0, label = str(gridsize_inverted[i]) +":")
            pl, =ax2.plot(data[0,:],data[1,:],col[i],linestyle ="solid", label = names[i]+ " $\u03C8_{s, bulk}$") # mean bulk soil pressure
            pl, =ax2.plot(data[0,:],data[2,:],col[i], linestyle ="dashed",label = names[i]+" $\u03C8_{s, eq}$") # ESWP
            #plt.plot(data[0,:],data[3,:],col[i],linestyle ="solid", label = "root collar pressure") # root collar pressure
            gridsize = ["1.0 cm", "multiscale"] #,  "1.0 cm Rhizo"]
            gridsize_inverted = gridsize[::-1]
            
            #plt.legend(ncol = 5, prop ={"size":9}, loc='lower center')

        except Exception as ex:
            print("Something went wrong with file "+f)    
            print(ex)

plt.plot(data_a_act[:,0], ESWP_reference, "k", linestyle="dashed", label="reference $\psi_{s, eq}$")
plt.plot(data_a_act[:,0], -head, "k", linestyle="dotted", label="reference $\psi_{s, bulk}$")

data_tcum2 = np.loadtxt("cumT_plots2/27-05-2021_CStudy_classic_sink_loam_dry_1,0cm",delimiter=';') # index 
tcum2 = data_tcum2[3,:]
t2 = data_tcum2[0,:] 

NR = (1 + (Alpha * abs(initialP))**N)**M
WC_i = Qr + (Qs - Qr)/NR

print("Initial water content: ", WC_i)

DomainWater = WC_i * 960
print("Amount of water in domain: ", DomainWater)

WC2 = (DomainWater - tcum2)/960
print(WC)
head2 = 1/Alpha*(((Qs-Qr)/(WC2 - Qr))**(1/M) - 1)**(1/N)

ax2.lines.pop(0)
ax2.plot(t2, -head2, color="b",linestyle ="solid", label = "1.0 cm RSS"+ " $\u03C8_{s, bulk}$") 
ax2.set_title("(B)")
ax2.set_xlabel("Time [days]")  
ax2.set_ylabel("Water potential [cm]")

handles, labels = plt.gca().get_legend_handles_labels()
order = [5,1,4,0,2,3]
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order], ncol=2, loc="lower center", fontsize=9 ) 
#plt.legend(loc="lower center", ncol=2, fontsize=9)



plt.xlim(0,x)
plt.ylim(-19000,0) # default is -19.000, 0!

plt.subplots_adjust(wspace=0.25)

#ax1.set_title('Tact Convergence Study Clay dry')plt.tight_layout()
plt.savefig(scenario+"_ESWP__rhizo_classic" + str(gridsize_inverted[i]) + "_same_color.jpg", dpi=300, bbox_inches = "tight")
plt.show()
            



'''
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
#gridsize = ["4.0 cm", "3.0 cm", "2.0 cm", "1.5 cm", "1.0 cm", "0.8 cm", "0.6 cm", "0.4 cm"] #,  "1.0 cm Rhizo"]
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
            #gridsize = ["4.0 cm", "3.0 cm", "2.0 cm", "1.5 cm", "1.0 cm", "0.8 cm", "0.6 cm", "0.4 cm"] #,  "1.0 cm Rhizo"]
            gridsize_inverted = gridsize[::-1]
            #plt.legend(ncol = 5, prop ={"size":9}, loc='lower center')
            plt.legend(bbox_to_anchor=(0,-0.2,1,0.2), loc="lower left",mode="expand", borderaxespad=0, ncol=4)

            plt.xlim(0,x)
            plt.ylim(-15000, 0)
            #ax1.set_title('Tact Convergence Study Clay dry')
            plt.tight_layout()
            plt.savefig(scenario+"_ESWP_" + str(gridsize_inverted[i]) + "_same_color.png", dpi=300, bbox_inches = "tight")
            #plt.show()
            


            plt.plot(data[0,:],data[1,:],"brown",linestyle ="solid", label = "bulk soil pressure") # mean bulk soil pressure
            plt.plot(data[0,:],data[2,:],"b", linestyle ="dashed", label = "ESWP") # ESWP
            plt.plot(data[0,:],data[3,:],"g",linestyle ="dotted", label = "root collar pressure") # root collar pressure
            plt.xlabel("Time [days]")  
            plt.ylabel("pressure head [cm]")
            #gridsize = ["4.0 cm", "3.0 cm", "2.0 cm", "1.5 cm", "1.0 cm", "0.8 cm", "0.6 cm", "0.4 cm"] #,  "1.0 cm Rhizo"]
            gridsize_inverted = gridsize[::-1]
            plt.legend(bbox_to_anchor=(0,-0.2,1,0.2), loc="lower left",mode="expand", borderaxespad=0, ncol=4)

            plt.xlim(0,x)
            plt.ylim(-15000, 0)
            #ax1.set_title('Tact Convergence Study Clay dry')
            plt.tight_layout()
            plt.savefig(scenario+"_ESWP_" + str(gridsize_inverted[i]) + "_different_color.png", dpi=300, bbox_inches = "tight")
            #plt.show()

        except Exception as ex:
            print("Something went wrong with file "+f)    
            print(ex)
'''
