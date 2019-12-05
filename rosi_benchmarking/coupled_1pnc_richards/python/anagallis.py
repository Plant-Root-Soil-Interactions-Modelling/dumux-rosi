''' Test '''

import os
import matplotlib.pyplot as plt
from vtk_tools import *
import van_genuchten as vg
import math

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/rosi_benchmarking/coupled_1pnc_richards")

# # run simulation
os.system("./coupled_1pnc_richards input/anagallis.input")

#      * 0 time [s], 1 actual transpiration [kg/s], 2 potential transpiration [kg/s], 3 maximal transpiration [kg/s],
#      * 4 collar pressure [Pa], 5 calculated actual transpiration [cm^3/day], 6 simtime [s], 7 hormone leaf mass [kg],
#      * 8 hormone collar flow rate [kg/s], 9 hormone root system mass [kg] , 10 hormone source rate [kg/s]
with open("anagallis_actual_transpiration.txt", 'r') as f:
    d = np.loadtxt(f, delimiter = ',')
c = 24 * 3600  # s / day

""" Plot transpiration """
plt.plot(d[:, 0] / c, 1000 * d[:, 2] * c, 'k')  # potential transpiration
plt.plot(d[:, 0] / c, 1000 * d[:, 1] * c, 'r-,')  # actual transpiration
plt.xlabel("time (days)")
plt.ylabel("transpiration (g/day)")
plt.legend(["potential", "actual"])
plt.title("Water transpiration")

""" Plot hormone rate and mass """
# fig, [ax1, ax2] = plt.subplots(1, 2)
#
# ax1.plot(d[:, 0] / c, 1000 * d[:, 8] * c, "r")
# ax1.plot(d[:, 0] / c, 1000 * d[:, 10] * c, "b")
# ax1.set_ylabel("mass rate (g/day)")
# ax1.set_xlabel("time (days)")
# ax1.set_title("Hormone production rate")
# ax1.legend(["leaf rate", "root system rate"])
#
# ax2.plot(d[:, 0] / c, 1000 * d[:, 7], "r")
# ax2.plot(d[:, 0] / c, 1000 * d[:, 9], "b")
# ax2.set_ylabel("mass (g)")
# ax2.set_xlabel("time (days)")
# ax2.set_title("Hormone mass")
# ax2.legend(["leaf mass", "root system mass"])
