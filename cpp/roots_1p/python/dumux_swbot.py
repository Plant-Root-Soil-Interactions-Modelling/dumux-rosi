#
# Wet bot scenario for 3 phenotypes in static soil (TODO Reference DGF missing)
#
""" reference dgf is missing """
import sys; sys.path.append("../../../../CPlantBox/src/python_modules/")

import os
import matplotlib.pyplot as plt
from vtk_tools import *
from math import *
import van_genuchten as vg
import threading
import time


class myThread(threading.Thread):

    def __init__(self, name):
        threading.Thread.__init__(self)
        self.name = name

    def run(self):
        print("Starting:" + self.name)
        os.system(self.name)


# Go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/cpp/roots_1p")

# Run dumux
t0 = time.time()
threads_ = []
threads_.append(myThread("./rootsystem input/swbot.input -RootSystem.Grid.File grids/RootSys/Reference/RootSys1.dgf"))
threads_.append(myThread("./rootsystem input/swbot.input -RootSystem.Grid.File grids/RootSys/Genotype_laterals/RootSys1.dgf -Problem.Name swbot_b"))
threads_.append(myThread("./rootsystem input/swbot.input -RootSystem.Grid.File grids/RootSys/Genotype_volume/RootSys1.dgf -Problem.Name swbot_c"))
for t in threads_:  # start threads
    t.start()
for t in threads_:  # and for all of them to finish
     t.join()
print("elapsed time is ", time.time() - t0)

# Read results
with open("swbot_actual_transpiration.txt", 'r') as f:
    d = np.loadtxt(f, delimiter = ',')
with open("swbot_b_actual_transpiration.txt", 'r') as f:
    d2 = np.loadtxt(f, delimiter = ',')
with open("swbot_c_actual_transpiration.txt", 'r') as f:
    d3 = np.loadtxt(f, delimiter = ',')

# Plot collar transpiration & pressure
fig, ax1 = plt.subplots()

c = (24 * 3600) / (.75 * .15)  #  -> mm / per day
t = d[:, 0] / (24 * 3600)
t2 = d2[:, 0] / (24 * 3600)
t3 = d3[:, 0] / (24 * 3600)

# 0 time, 1 actual transpiration, 2 potential transpiration, 3 maximal transpiration, 4 collar pressure, 5 calculated actual transpiration
ax1.plot(t, d[:, 2] * c, 'k')  # potential transpiration
ax1.plot(t, d[:, 1] * c, 'r-,')  # reference, actual transpiration
ax1.plot(t, d[:, 3] * c, 'r:')  # reference, maximal transpiration
ax1.plot(t2, d2[:, 1] * c, 'g-,')  # lateral
ax1.plot(t2, d2[:, 3] * c, 'g:')  # lateral
ax1.plot(t3, d3[:, 1] * c, 'b-,')  # volume
ax1.plot(t3, d3[:, 3] * c, 'b:')  # volume
ax1.legend(['Pot trans', 'actual trans P1', 'max trans P1', 'actual trans P2', 'max trans P2', 'actual trans P3', 'max trans P3'], loc = 'upper left')
ax1.axis((0, t[-1], 0, 12))
ax1.set_xlabel("Time $[d]$")
ax1.set_ylabel("Transpiration rate $[mm \ d^{-1}]$")

# ax2 = ax1.twinx()
# ax2.plot(t, d[:, 5], 'r--')
# ax2.plot(t2, d2[:, 5], 'g--')
# # ax2.plot(t, d3[:, 5], 'b--')
# ax2.legend(['pressure P1', 'pressure P2', 'pressure P3'], loc = 'upper right')
# ax2.set_ylabel("Pressure at root collar (Pa)")

# Some text output
i = np.argmax(d[:, 1] == d[:, 3])
print("stress after ", d[i, 0] / (24 * 3600), "days")
i = np.argmax(d2[:, 1] == d2[:, 3])
print("stress after ", d[i, 0] / (24 * 3600), "days")
i = np.argmax(d3[:, 1] == d3[:, 3])
print("stress after ", d[i, 0] / (24 * 3600), "days")

print("transpiration during stress", d[-1, 1] * c, "kg/s at", d[-1, 4], "Pa, p-crit ", d[-1, 5], "Pa, max trans", d[-1, 3])
print("transpiration during stress", d2[-1, 1] * c, "kg/s at", d2[-1, 4], "Pa, p-crit ", d2[-1, 5], "Pa, max trans", d2[-1, 3])
print("transpiration during stress", d3[-1, 1] * c, "kg/s at", d3[-1, 4], "Pa, p-crit ", d3[-1, 5], "Pa, max trans", d3[-1, 3])

plt.show()

