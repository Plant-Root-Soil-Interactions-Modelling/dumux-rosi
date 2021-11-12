''' Soybean stomatal coupled model '''
import sys; sys.path.append("../../../../CPlantBox/src/python_modules/")

import os
import matplotlib.pyplot as plt
from vtk_tools import *
import van_genuchten as vg
import math

name = "soybean_Conly_2010"  # this name should be unique
suffix = ""

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/cpp/coupled_1pnc_richards")

# run simulation
os.system("./coupled_periodic_1pnc_richards input/" + name + ".input")

# move results to folder 'name'
if not os.path.exists("results_" + name + suffix):
    os.mkdir("results_" + name + suffix)
os.system("mv " + name + "* " + "results_" + name + suffix + "/")
os.system("cp input/" + name + ".input " + "results_" + name + suffix + "/")

#      * 0 time [s], 1 actual transpiration [kg/s], 2 potential transpiration [kg/s], 3 maximal transpiration [kg/s],
#      * 4 collar pressure [Pa], 5 calculated actual transpiration [cm^3/day], 6 simtime [s], 7 hormone leaf mass [kg],
#      * 8 hormone collar flow rate [kg/s], 9 hormone root system mass [kg] , 10 hormone source rate [kg/s], 

with open("results_" + name + suffix + "/" + name + "_actual_transpiration.txt", 'r') as f:
    d = np.loadtxt(f, delimiter = ',')

c = 24 * 3600  # s / day
t = d[:, 0] / c  # [s] -> [day]

""" Plot hormone rate and mass """
fig, ax1 = plt.subplots()

ax1.plot(d[:, 0] / c, 1000 * d[:, 8] * c, "r")
ax1.tick_params(axis = 'y', labelcolor = "r")
ax1.set_ylabel("Leaves [g/day]", color = "r")
ax1b = ax1.twinx()
ax1b.plot(d[:, 0] / c, 1000 * d[:, 10] * c, "b")
ax1b.tick_params(axis = 'y', labelcolor = "b")
ax1b.set_ylabel("Root system [g/day]", color = "b")
ax1.set_xlabel("Time [days]")
ax1.set_title("Hormone production rate")
plt.savefig("results_" + name + suffix + "/" + name + suffix + '_productionrate.pdf', dpi=300, bbox_inches='tight')

fig, ax2 = plt.subplots()

ax2.plot(d[:, 0] / c, 1000 * d[:, 7], "r")
ax2.tick_params(axis = 'y', labelcolor = "r")
ax2.set_ylabel("Leaves [g]", color = "r")
ax2b = ax2.twinx()
ax2b.plot(d[:, 0] / c, 1000 * d[:, 9], "b")
ax2b.tick_params(axis = 'y', labelcolor = "b")
ax2b.set_ylabel("Root system [g]", color = "b")
ax2.set_xlabel("Time [days]")
ax2.set_title("Hormone mass")
plt.savefig("results_" + name + suffix + "/" + name + suffix + '_hormoneMass.pdf', dpi=300, bbox_inches='tight')

""" Plot transpiration """
fig, ax3 = plt.subplots()
ax3.plot(t, 1000 * d[:, 2] * c, 'k')  # potential transpiration
ax3.plot(t, 1000 * d[:, 1] * c, 'r-,')  # actual transpiration
ax3b = ax3.twinx()
ctrans = np.cumsum(np.multiply(1000 * d[1:, 1] * c, (t[1:] - t[:-1])))
ax3b.plot(t[1:], ctrans, 'c--', color = 'blue')  # cumulative transpiration (neumann)
ax3b.tick_params(axis= 'y', labelcolor = 'b')
ax3.set_xlabel("Time [days]")
ax3.set_ylabel("Transpiration [mL/day]")
ax3b.set_ylabel("Cumulative transpiration $[mL]$", color = "b")
ax3.legend(["potential", "actual", "cumulative"], loc = 'upper left')
plt.savefig("results_" + name + suffix + "/" + name + suffix + '_transpiration.pdf', dpi=300, bbox_inches='tight')

""" Plot root collar pressure """
fig, ax4 = plt.subplots()
ax4.plot(d[:, 0] / c, (d[:, 4] - 1.e5) * 100. / 1.e3 / 9.81, 'r-')  # root collar pressure head (convert from Pa to Head)  
ax4.set_xlabel("Time [days]")
ax4.set_ylabel("Pressure at root collar [cm]")
plt.savefig("results_" + name + suffix + "/" + name + suffix + '_collarPressure.pdf', dpi=300, bbox_inches='tight')

plt.show()
