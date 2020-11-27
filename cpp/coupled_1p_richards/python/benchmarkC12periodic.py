''' Run Benchmakrk C12periodic (like benchmarkC12.py but with periodic boundary conditions)'''
import sys; sys.path.append("../../../../CPlantBox/src/python_modules/")

import os
import matplotlib.pyplot as plt
from vtk_tools import *
import van_genuchten as vg
import math

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/cpp/coupled_1p_richards")

# run simulation
os.system("./coupled_periodic input/benchmarkC12periodic.input")

# 0 time [s], 1 actual transpiration [kg/s], 2 potential transpiration [kg/s], 3 maximal transpiration [kg/s],
# 4 collar pressure [Pa], 5 calculated actual transpiration, 6 time [s]
with open("benchmarkC12periodic_actual_transpiration.txt", 'r') as f:  # benchmarkC12c_actual_transpiration. or benchmarkC12bc_actual_transpiration
    d = np.loadtxt(f, delimiter = ',')

print()
c = 24 * 3600  #  [kg/s] -> [kg/per day]
print("potential", d[-1, 2] * c)
print("actual", d[-1, 1] * c)
print("actual", d[-1, 5] / 1000)  # Strange behaviour of simplistically calculated radial flows

# Plot collar transpiration & pressure
fig, ax1 = plt.subplots()

c = 1000 * 24 * 3600  #  [kg/s] -> [cm3/per day]
t = d[:, 0] / (24 * 3600)  # [s] -> [day]

# 0 time, 1 actual transpiration, 2 potential transpiration, 3 maximal transpiration, 4 collar pressure, 5 calculated actual transpiration
ax1.plot(t, d[:, 2] * c, 'k')  # potential transpiration
ax1.plot(t, d[:, 1] * c, 'g')  # actual transpiration (neumann)

ax2 = ax1.twinx()
ctrans = np.cumsum(np.multiply(d[1:, 1] * c, (t[1:] - t[:-1])))
ax2.plot(t[1:], ctrans, 'c--')  # cumulative transpiration (neumann)

ax1.legend(['Potential', 'Actual', 'Cumulative'], loc = 'upper left')
ax1.set_xlabel("Time $[d]$")
ax1.set_ylabel("Transpiration rate $[cm^3 \ d^{-1}]$")
ax2.set_ylabel("Cumulative transpiration $[cm^3]$")

plt.show()

# trans = interpolate.interp1d(t, d[:, 1] * c)
# print(t.shape)
# ctrans = np.zeros(t.shape)
# ctrans[0] = 0
# for i, t_ in enumerate(t[1:]):
#     ctrans[i] = integrate.quad(trans, 0, t_)[0]
