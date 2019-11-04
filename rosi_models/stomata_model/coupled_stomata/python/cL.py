''' Test '''

import os
import matplotlib.pyplot as plt
from vtk_tools import *
import van_genuchten as vg
import math

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../../build-cmake/rosi_models/stomata_model/coupled_stomata")

with open("cL.txt", 'r') as f:  # benchmarkC12c_actual_transpiration. or benchmarkC12bc_actual_transpiration
    d = np.loadtxt(f, delimiter = ',')

# Plot chemical concentration
fig, ax1 = plt.subplots()

t = d[:, 0] / (24 * 3600)  # [s] -> [day]

# 0 time, 1 actual transpiration, 2 potential transpiration, 3 maximal transpiration, 4 collar pressure, 5 calculated actual transpiration, 6 time, 7 chemical concentration
ax1.plot(t, d[:, 1], 'k')  # potential transpiration

ax1.legend(['chemical concentration'], loc = 'upper left')
#ax1.axis((0, t[-1], 0, 13))
ax1.set_xlabel("Time $[d]$")
ax1.set_ylabel("chemical concentration $[cm^3 \ d^{-1}]$")

#plt.savefig('../results/Transpitation_2000.png')
plt.show()

# trans = interpolate.interp1d(t, d[:, 1] * c)
# print(t.shape)
# ctrans = np.zeros(t.shape)
# ctrans[0] = 0
# for i, t_ in enumerate(t[1:]):
#     ctrans[i] = integrate.quad(trans, 0, t_)[0]
