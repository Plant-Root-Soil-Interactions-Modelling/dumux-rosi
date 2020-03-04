''' singleroot '''

import os
import matplotlib.pyplot as plt
from vtk_tools import *
import van_genuchten as vg
import math

name = "singleroot"  # this name should be unique

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/rosi_benchmarking/coupled_1p_richards")

# # run simulation
os.system("./coupled_periodic input/singleroot.input")

# move results to folder 'name'
if not os.path.exists("results_" + name):
    os.mkdir("results_" + name)
os.system("mv " + name + "* " + "results_" + name + "/")
os.system("cp input/" + name + ".input " + "results_" + name + "/")

# 0 time [s], 1 actual transpiration [kg/s], 2 potential transpiration [kg/s], 3 maximal transpiration [kg/s],
# 4 collar pressure [Pa], 5 calculated actual transpiration, 6 time [s]
with open("results_" + name + "/" + name + "_actual_transpiration.txt", 'r') as f:  # benchmarkC12c_actual_transpiration. or benchmarkC12bc_actual_transpiration
    d = np.loadtxt(f, delimiter = ',')
c = 24 * 3600  #  [kg/s] -> [kg/per day]
t = d[:, 0] / c  # [s] -> [day]

print("potential", d[-1, 2] * c)
print("actual", d[-1, 1] * c)
print("actual", d[-1, 5] / 1000)  # Strange behaviour of simplistically calculated radial flows

# Plot collar transpiration & pressure
fig, ax1 = plt.subplots()

# 0 time, 1 actual transpiration, 2 potential transpiration, 3 maximal transpiration, 4 collar pressure, 5 calculated actual transpiration
ax1.plot(t, 1000 * d[:, 2] * c, 'k')  # potential transpiration
ax1.plot(t, 1000 * d[:, 1] * c, 'r-')  # actual transpiration (neumann)

ax2 = ax1.twinx()
ctrans = np.cumsum(np.multiply(1000 * d[1:, 1] * c, (t[1:] - t[:-1])))
ax2.plot(t[1:], ctrans, 'c--', color = 'blue')  # cumulative transpiration (neumann)
ax2.tick_params(axis= 'y', labelcolor = 'b')
ax1.set_xlabel("time (days)")
ax1.set_ylabel("transpiration rate (g/day)")
ax2.set_ylabel("Cumulative transpiration $[g]$", color = "b")
ax1.legend(['Potential', 'Actual', 'Cumulative'], loc = 'upper left')
plt.savefig("results_" + name + ".pdf", dpi=300)
plt.show()

