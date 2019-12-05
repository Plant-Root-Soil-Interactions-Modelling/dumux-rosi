""" "Test root system part """
import os
import matplotlib.pyplot as plt
from vtk_tools import *
from math import *
import van_genuchten as vg
import threading
import time

# Go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/rosi_benchmarking/roots_1pnc")

# run dumux
# os.system("./rootsystem_1pnc input/anagallis_stomata.input")

#      * 0 time [s], 1 actual transpiration [kg/s], 2 potential transpiration [kg/s], 3 maximal transpiration [kg/s],
#      * 4 collar pressure [Pa], 5 calculated actual transpiration [cm^3/day], 6 simtime [s], 7 hormone leaf mass [kg],
#      * 8 hormone collar flow rate [kg/s], 9 hormone root system mass [kg] , 10 hormone source rate [kg/s]
with open("anagallis_roots_actual_transpiration.txt", 'r') as f:
    d = np.loadtxt(f, delimiter = ',')
c = 24 * 3600  # s / day

# """ Plot transpiration """
# plt.plot(d[:, 0] / c, 1000 * d[:, 2] * c, 'k')  # potential transpiration
# plt.plot(d[:, 0] / c, 1000 * d[:, 1] * c, 'r-,')  # actual transpiration
# plt.xlabel("time (days)")
# plt.ylabel("transpiration (g/day)")
# plt.legend(["potential", "actual"])
# plt.title("Water transpiration")

""" Plot hormone rate and mass """
fig, [ax1, ax2] = plt.subplots(1, 2)

ax1.plot(d[:, 0] / c, 1000 * d[:, 8] * c, "r")
ax1.tick_params(axis = 'y', labelcolor = "r")
ax1.set_ylabel("leaf mass rate (g/day)", color = "r")
ax1b = ax1.twinx()
ax1b.plot(d[:, 0] / c, 1000 * d[:, 10] * c, "b")
ax1b.tick_params(axis = 'y', labelcolor = "b")
ax1b.set_ylabel("root system production rate (g/day)", color = "b")
ax1.set_xlabel("time (days)")
ax1.set_title("Hormone production rate")

ax2.plot(d[:, 0] / c, 1000 * d[:, 7], "r")
ax2.tick_params(axis = 'y', labelcolor = "r")
ax2.set_ylabel("leaf mass (g)", color = "r")
ax2b = ax2.twinx()
ax2b.plot(d[:, 0] / c, 1000 * d[:, 9], "b")
ax2b.tick_params(axis = 'y', labelcolor = "b")
ax2b.set_ylabel("root system mass (g)", color = "b")
ax2.set_xlabel("time (days)")
ax2.set_title("Hormone mass")

plt.show()
