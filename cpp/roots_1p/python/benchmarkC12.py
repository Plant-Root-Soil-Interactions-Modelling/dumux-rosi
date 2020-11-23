""" "Benchmark C12 root system part """
import sys; sys.path.append("../../../../CPlantBox/src/python_modules/")

import os
import matplotlib.pyplot as plt
from vtk_tools import *
import van_genuchten as vg

# Go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/cpp/roots_1p")

# run dumux
os.system("./rootsystem input/benchmarkC12.input")

# 0 time [s], 1 actual transpiration [kg/s], 2 potential transpiration [kg/s], 3 maximal transpiration [kg/s],
# 4 collar pressure [Pa], 5 calculated actual transpiration, 6 time [s]
with open("benchmarkC12_actual_transpiration.txt", 'r') as f:
    d = np.loadtxt(f, delimiter = ',')

print()
c = 24 * 3600  #  [kg/s] -> [kg/per day]
print("potential", d[-1, 2] * c)
print("actual", d[-1, 1] * c)
print("actual", d[-1, 5] / 1000)

# Plot collar transpiration & pressure
fig, ax1 = plt.subplots()

c = 24 * 3600  #  [kg/s] -> [kg/per day]
t = d[:, 6] / (24 * 3600)  # [s] -> [day]

# 0 time, 1 actual transpiration, 2 potential transpiration, 3 maximal transpiration, 4 collar pressure, 5 calculated actual transpiration
ax1.plot(t, d[:, 2] * c, 'k')  # potential transpiration
ax1.plot(t, d[:, 5] / 1000, 'r-,')  # actual transpiration (calculated)
ax1.plot(t, d[:, 1] * c, 'g:')  # actual transpiration (neumann)

ax1.legend(['Potential', 'Actual', 'Actual'], loc = 'upper left')
ax1.axis((0, t[-1], 0, 0.013))
ax1.set_xlabel("Time $[d]$")
ax1.set_ylabel("Transpiration rate $[kg \ d^{-1}]$")

plt.show()
