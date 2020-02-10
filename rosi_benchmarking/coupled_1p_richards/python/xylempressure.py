''' singleroot '''

import os
import matplotlib.pyplot as plt
from vtk_tools import *
import van_genuchten as vg
import math
import numpy as np

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/rosi_benchmarking/coupled_1pnc_richards")

#      * 0 time [s], 1 actual transpiration [kg/s], 2 potential transpiration [kg/s], 3 maximal transpiration [kg/s],
#      * 4 collar pressure [Pa], 5 calculated actual transpiration [cm^3/day], 6 simtime [s], 7 hormone leaf mass [kg],
#      * 8 hormone collar flow rate [kg/s], 9 hormone root system mass [kg] , 10 hormone source rate [kg/s]
with open("singleroot_actual_transpiration.txt", 'r') as f:
    d = np.loadtxt(f, delimiter = ',')

""" Plot Tact/Tpot vs xylem pressure"""
p_, z_ = read3D_vtp_data("singleroot-00058.vtp")
h_ = vg.pa2head(p_)
#plt.plot(h_, d[:, 1]/d[:, 2], "r+")
plt.plot(h_, z_[:,2], "r+")  # cell data
plt.set_ylabel("$alpha$")
plt.set_xlabel("Xylem pressure (cm)")

plt.show()
