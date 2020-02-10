''' singleroot '''

import os
import matplotlib.pyplot as plt
from vtk_tools import *
import van_genuchten as vg
import math

#      * 0 time [s], 1 actual transpiration [kg/s], 2 potential transpiration [kg/s], 3 maximal transpiration [kg/s],
#      * 4 collar pressure [Pa], 5 calculated actual transpiration [cm^3/day], 6 simtime [s], 7 hormone leaf mass [kg],
#      * 8 hormone collar flow rate [kg/s], 9 hormone root system mass [kg] , 10 hormone source rate [kg/s]
with open("singleroot_actual_transpiration.txt", 'r') as f:
    d = np.loadtxt(f, delimiter = ',')
c = 24 * 3600  # s / day
t = d[:, 0] / c  # [s] -> [day]

""" Plot Tact/Tpot """
fig, ax1 = plt.subplots()
ax1.plot(d[:, 0] / c, 1000 * (d[:, 1]/d[:, 2]) * c, 'r-')  # actual transpiration/potential transpiration
ax1.set_xlabel("time (days)")
ax1.set_ylabel("T$_{act}/T_{pot}$")
ax1.set_title("Stress (T$_{act}/T_{pot}$)")

plt.show()
