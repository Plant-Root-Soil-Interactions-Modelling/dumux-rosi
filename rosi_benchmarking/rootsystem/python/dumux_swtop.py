#
# Wet top scenario for 3 phenotypes in static soil
#

import os
import matplotlib.pyplot as plt
from vtk_tools import *
from math import *
import van_genuchten as vg
import time

g = 9.81  # gravitational acceleration (m/s^2)
rho = 1.e3  # density of water, (kg/m^3)
ref = 1.e5


def toHead(pa):  # Pascal (kg/ (m s^2)) to cm pressure head
    return (pa - ref) * 100 / rho / g


# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/rosi_benchmarking/rootsystem")

print("wet", toHead(-10000))
print("dry", toHead(-300000))
print("a week ", 7 * 24 * 3600)

trans = 5.33 * .75 * .15 / 86400
maxtrans = 2 * trans
print("daily rate ", 5.33, "mm/day = ", trans, " kg/s, maximum ", maxtrans)  #

# run dumux
t = time.time()
os.system("./rootsystem input/swtop.input -RootSystem.Grid.File grids/RootSys/Reference/RootSys1.dgf")  # mpirun -n 8
os.system("./rootsystem input/swtop.input -RootSystem.Grid.File grids/RootSys/Genotype_laterals/RootSys1.dgf -Problem.Name swtop_b")  # mpirun -n 8
os.system("./rootsystem input/swtop.input -RootSystem.Grid.File grids/RootSys/Genotype_volume/RootSys1.dgf -Problem.Name swtop_c")  # mpirun -n 8
elapsed = time.time() - t
print("elapsed time is ", elapsed)

with open("swtop_actual_transpiration.txt", 'r') as f:
    d = np.loadtxt(f, delimiter = ',')
with open("swtop_b_actual_transpiration.txt", 'r') as f:
    d2 = np.loadtxt(f, delimiter = ',')
with open("swtop_c_actual_transpiration.txt", 'r') as f:
    d3 = np.loadtxt(f, delimiter = ',')
# Format of txt file:
# time_, lastActualTrans_, lastTrans_, lastMaxTrans_, (sol[0] - pRef_) * 100 / rho_ / g_, trans
# 0    , 1               , 2         , 3            , 4                                 , 5
#

# days = 1
# t_ = np.linspace(0, days, 6 * 24 * days)
# y_ = np.sin(t_ * 2.*pi - 0.5 * pi) * trans + trans
# ax1.plot(t_, y_ * (24 * 3600), 'k')

# Plot collar transpiration & pressure
fig, ax1 = plt.subplots()

ax1.plot(d[:, 0] / (24 * 3600), d[:, 2] * (24 * 3600), 'k')  # potential transpiration
ax1.plot(d[:, 0] / (24 * 3600), d[:, 1] * (24 * 3600), 'r')  # reference, actual transpiration
ax1.plot(d[:, 0] / (24 * 3600), d[:, 3] * (24 * 3600), 'r:')
ax1.plot(d2[:, 0] / (24 * 3600), d2[:, 1] * (24 * 3600), 'g')  # lateral
ax1.plot(d2[:, 0] / (24 * 3600), d2[:, 3] * (24 * 3600), 'g:')  # lateral
ax1.plot(d3[:, 0] / (24 * 3600), d3[:, 1] * (24 * 3600), 'b')  # volume
ax1.plot(d3[:, 0] / (24 * 3600), d3[:, 1] * (24 * 3600), 'b:')  # volume
ax1.legend(['Pot trans', 'actual trans P1', 'max trans P1', 'actual trans P2', 'max trans P2', 'actual trans P3', 'max trans P3'], loc = 'upper left')
ax1.axis((0, days, 0, 1.3))
ax1.set_xlabel("Time $[d]$")
ax1.set_ylabel("Transpiration rate $[kg d^{-1}]$")

ax2 = ax1.twinx()
ax2.plot(d[:, 0] / (24 * 3600), d[:, 4], 'r--')
ax2.plot(d2[:, 0] / (24 * 3600), d2[:, 4], 'g--')
ax2.plot(d3[:, 0] / (24 * 3600), d3[:, 4], 'b--')
ax2.legend(['pressure P1', 'pressure P2', 'pressure P3'], loc = 'upper right')

ax2.set_ylabel("Pressure at root collar (cm)")
print(np.min(d[:, 4]))

plt.show()

i = np.argmax(d[:, 1] == d[:, 3])
print("stress after ", d[i, 0] / (24 * 3600), "days")
i = np.argmax(d2[:, 1] == d2[:, 3])
print("stress after ", d[i, 0] / (24 * 3600), "days")
i = np.argmax(d3[:, 1] == d3[:, 3])
print("stress after ", d[i, 0] / (24 * 3600), "days")

# for 1 hour time step
# stress after  0.396412037037 days
# stress after  0.146412037037 days
# stress after  0.146412037037 days

# # plot
# p_, z_ = read3D_vtp_data("swtop-00001.vtp", False)
# h_ = vg.pa2head(p_)
# plt.plot(p_, z_[:, 2], "r+")  # cell data
# plt.ylabel("Depth (m)")
# plt.xlabel("Xylem pressure (cm)")
# plt.show()

