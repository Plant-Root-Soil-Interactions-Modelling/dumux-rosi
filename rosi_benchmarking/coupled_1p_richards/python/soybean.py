''' Soybean model '''

import os
import matplotlib.pyplot as plt
from vtk_tools import *
import van_genuchten as vg
import math

name = "soybean_Honly"  # this name should be unique
suffix = ""

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/rosi_benchmarking/coupled_1p_richards")

# run simulation
os.system("./coupled_periodic input/" + name + ".input -RootSystem.Grid.InitialT 1")

# move results to folder 'name'
#if not os.path.exists("results_" + name + suffix):
#    os.mkdir("results_" + name + suffix)
#os.system("mv " + name + "* " + "results_" + name + suffix + "/")
#os.system("cp input/" + name + ".input " + "results_" + name + suffix + "/")

# 0 time [s], 1 actual transpiration [kg/s], 2 potential transpiration [kg/s], 3 maximal transpiration [kg/s],
# 4 collar pressure [Pa], 5 calculated actual transpiration, 6 time [s]
with open("results_" + name + suffix + "/" + name + "_actual_transpiration.txt", 'r') as f:
    d = np.loadtxt(f, delimiter = ',')
c = 24 * 3600  # s / day
t = d[:, 0] / c  # [s] -> [day]

""" Plot transpiration """
fig, ax1 = plt.subplots()
ax1.plot(t, 1000 * d[:, 2] * c, 'k')  # potential transpiration
ax1.plot(t, 1000 * d[:, 1] * c, 'r-,')  # actual transpiration

ax1b = ax1.twinx()
ctrans = np.cumsum(np.multiply(1000 * d[1:, 1] * c, (t[1:] - t[:-1])))
ax1b.plot(t[1:], ctrans, 'c--', color = 'blue')  # cumulative transpiration (neumann)
ax1b.tick_params(axis= 'y', labelcolor = 'b')
ax1.set_xlabel("Time [days]")
ax1.set_ylabel("Transpiration [mL/day]")
ax1b.set_ylabel("Cumulative transpiration $[mL]$", color = "b")
ax1.legend(["potential", "actual", "cumulative"], loc = 'upper left')
plt.savefig("results_" + name + suffix + "/" + name + suffix + '_transpiration.pdf', dpi=300)

""" Plot root collar pressure """
fig, ax2 = plt.subplots()
ax2.plot(t, (d[:, 4] - 1.e5) * 100. / 1.e3 / 9.81, 'r-')  # root collar pressure head (convert from Pa to Head)  
ax2.set_xlabel("Time [days]")
ax2.set_ylabel("Pressure at root collar [cm]")
plt.savefig("results_" + name + suffix + "/" + name + suffix + '_collarPressure.pdf', dpi=300)

""" Equivalent soil water potential """
#fig, ax3 = plt.subplots()
#psi_eq = (1000 * d[:, 1] * c) + ((d[:, 4] - 1.e5) * 100. / 1.e3 / 9.81)
#ax3.plot(t, psi_eq)
#ax1.set_xlabel("Time [days]")
#ax1.set_ylabel("Equivalent soil water potential [cm]")
#s_, p_, z_ = read3D_vtp_data("results_" + name + suffix + "/" + "soybean_Honly-00480.vtu", 1)
#h_ = vg.pa2head(p_)
#plt.plot(h_, z_ * 100, "r+")
#plt.xlabel("Soil pressure (cm)")
#plt.ylabel("Depth (cm)")

plt.show()