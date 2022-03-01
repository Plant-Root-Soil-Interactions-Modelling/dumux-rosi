
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from math import *
import csv
from pathlib import Path

fig, ax1 = plt.subplots(figsize=(8,7)) 
col = ["r", "tab:pink","b", "c", "g"]

x = 7   # time for outputs
scenario = "clay_dry"

# ESWP
l=[]
for dirname, dirnames, filenames in os.walk('ESWP_SI/.'):
    filenames.sort()
    for i,f in enumerate(filenames):
        try:
            l.append(Path(f).stem)
            data = np.loadtxt("ESWP_SI/"+f,delimiter=';')
            names = ["1.0 cm multiscale", "0.4 cm RSS","1.0 cm RSS","2.0 cm RSS","4.0 cm RSS"]
            #pl, = ax1.plot(data[0,:],data[1,:],col[i]) # TPot
            #pl, = ax1.plot(data[0,:],data[2,:],col[i]) # ESWP
            plt.plot(data[0,:],data[2,:],col[i], linestyle="dashed", label = names[i]+" $\u03C8_{s, eq}$") # ESWP
            #plt.plot(data[0,:],data[1,:],col[i],linestyle ="solid", label = names[i]+ " $\u03C8_{s, bulk}$") # mean bulk soil pressure

            #pl, = ax1.plot(data[0,:],data[3,:],col[i]) # root collar pressure
            
        except Exception as ex:
            print("Something went wrong with file "+f)    
            print(ex)
#plt.xlabel("Time [days]")  
#plt.ylabel("Water potential [cm]")
#plt.legend( loc="best", fontsize=9) 
#plt.xlim(0,x)
#plt.ylim(-15000, 0)
#plt.tight_layout()
#plt.savefig(scenario+"_ESWP_all_SI.png", dpi=300, bbox_inches = "tight")
#print(l)
#plt.show()

l=[]
for dirname, dirnames, filenames in os.walk('ESWP_SI_cumT/.'):
    filenames.sort()
    for i,f in enumerate(filenames):
        try:
            l.append(Path(f).stem)
            data = np.loadtxt("ESWP_SI_cumT/"+f,delimiter=';')
            names = ["1.0 cm multiscale", "0.4 cm RSS","1.0 cm RSS","2.0 cm RSS","4.0 cm RSS"]
            t = data[0,:]
            Alpha = 0.01
            N = 1.1
            Qr = 0.1
            Qs = 0.4
            M = 1-1/N

            #clay = [0.1, 0.4, 0.01, 1.1, 10]

            initialP = -652.3
            NR = (1 + (Alpha * abs(initialP))**N)**M
            WC_i = Qr + (Qs - Qr)/NR

            print("Initial water content: ", WC_i)
            DomainWater = WC_i * 960
            print("Amount of water in domain: ", DomainWater)
            WC = (DomainWater - data[3,:])/960
            #print(WC)

            head = 1/Alpha*(((Qs-Qr)/(WC - Qr))**(1/M) - 1)**(1/N)
            plt.plot(t, -head,col[i],linestyle ="solid", label = names[i]+ " $\u03C8_{s, bulk}$") # mean bulk soil pressure

            #pl, = ax1.plot(data[0,:],data[3,:],col[i]) # root collar pressure
            
        except Exception as ex:
            print("Something went wrong with file "+f)    
            print(ex)
plt.xlabel("Time [days]")  
plt.ylabel("Water potential [cm]")
plt.legend( loc="best", fontsize=9) 
plt.xlim(0,x)
plt.ylim(-15000, 0)
plt.tight_layout()
plt.savefig(scenario+"_ESWP_all_SI.jpg", dpi=300, bbox_inches = "tight")
plt.show()

