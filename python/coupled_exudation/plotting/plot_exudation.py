"""
transpiration plot (one column, number of rows as number of filenames)
"""
import sys; sys.path.append("../../../CPlantBox/");  sys.path.append("../../../CPlantBox/src/python_modules")

import numpy as np
import matplotlib.pyplot as plt

from xylem_flux import sinusoidal2
import evapotranspiration as evap

Kc_maize = 1.2
Kc_soybean = 1.15

#
name = "maize"
str_ = "_cyl3"

fname = "carbon_" + name + str_
path = "results/"

#SMALL_SIZE = 16
#MEDIUM_SIZE = 16
#BIGGER_SIZE = 16
#plt.rc('font', size = SMALL_SIZE)  # controls default text sizes
#plt.rc('axes', titlesize = SMALL_SIZE)  # fontsize of the axes title
#plt.rc('axes', labelsize = MEDIUM_SIZE)  # fontsize of the x and y labels
#plt.rc('xtick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
#plt.rc('ytick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
#plt.rc('legend', fontsize = SMALL_SIZE)  # legend fontsize
#plt.rc('figure', titlesize = BIGGER_SIZE)  # fontsize of the figure title

""" exudation plot """

# load data
print(fname) 
data = np.load(path + fname + ".npy") 

fig, ax = plt.subplots(1, 1, figsize = (18, 8))

print(np.shape(data)) 
t = data[0]
c = data[1]
lns1 = ax.plot(t, c, 'r', label = "Carbon exudation [g/day]")
ax.set_ylabel("Carbon exudation [g/day]")
dt = np.diff(t)
cuc = np.cumsum(np.multiply(c[:-1], dt))
print(str_, "Cumulative carbon exudation", cuc[-1] , "g", c[0])

ax2 = ax.twinx()
ax2.set_ylabel("Cumulative carbon exudation [g]")
lns2 = ax2.plot(t[1:],  cuc, 'r--', label = "Cumulative carbon exudation [g]")  

# added these three lines
lns = lns1+lns2
labs = [l.get_label() for l in lns]
ax.legend(lns, labs, loc=0)
    
ax.set_xlabel("Time [day]")

#plt.tight_layout()
plt.show()
