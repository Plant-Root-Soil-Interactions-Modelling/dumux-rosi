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


def toPa(ph):
    return ref + ph / 100. * rho * g;


print("wet", toHead(-10000))
print("dry", toHead(-300000))
print("a week ", 7 * 24 * 3600)
trans = 5.33 * .75 * .15 / 86400
maxtrans = 2 * trans
print("daily rate ", 5.33, "mm/day = ", trans, " kg/s, maximum ", maxtrans)  #
print("Critical collar pressure = ", toPa(-1.e4))

# Go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/rosi_benchmarking/rootsystem")

# Run dumux
t = time.time()
os.system("./rootsystem input/swtop.input -RootSystem.Grid.File grids/RootSys/Reference/RootSys1.dgf")  # mpirun -n 8
# os.system("./rootsystem input/swtop.input -RootSystem.Grid.File grids/RootSys/Genotype_laterals/RootSys1.dgf -Problem.Name swtop_b")  # mpirun -n 8
os.system("./rootsystem input/swtop.input -RootSystem.Grid.File grids/RootSys/Genotype_volume/RootSys1.dgf -Problem.Name swtop_c")  # mpirun -n 8
print("elapsed time is ", time.time() - t)

with open("swtop_actual_transpiration.txt", 'r') as f:
    d = np.loadtxt(f, delimiter = ',')
with open("swtop_b_actual_transpiration.txt", 'r') as f:
    d2 = np.loadtxt(f, delimiter = ',')
with open("swtop_c_actual_transpiration.txt", 'r') as f:
    d3 = np.loadtxt(f, delimiter = ',')
# Format of txt file:
# time_, lastActualTrans_, lastTrans_, lastMaxTrans_, (p - pRef_) * 100 / rho_ / g_, dp, trans
# 0    , 1               , 2         , 3            , 4                            , 5 , 6
#

# Plot collar transpiration & pressure
fig, ax1 = plt.subplots()

ax1.plot(d[:, 0] / (24 * 3600), d[:, 2] * (24 * 3600), 'k')  # potential transpiration
ax1.plot(d[:, 0] / (24 * 3600), d[:, 1] * (24 * 3600), 'r')  # reference, actual transpiration
ax1.plot(d[:, 0] / (24 * 3600), d[:, 3] * (24 * 3600), 'r:')  # reference, maximal transpiration
ax1.plot(d2[:, 0] / (24 * 3600), d2[:, 1] * (24 * 3600), 'g')  # lateral
ax1.plot(d2[:, 0] / (24 * 3600), d2[:, 3] * (24 * 3600), 'g:')  # lateral
ax1.plot(d3[:, 0] / (24 * 3600), d3[:, 1] * (24 * 3600), 'b')  # volume
ax1.plot(d3[:, 0] / (24 * 3600), d3[:, 1] * (24 * 3600), 'b:')  # volume
ax1.legend(['Pot trans', 'actual trans P1', 'max trans P1', 'actual trans P2', 'max trans P2', 'actual trans P3', 'max trans P3'], loc = 'upper left')
ax1.axis((0, d[-1, 0] / (24 * 3600), 0, 1.3))
ax1.set_xlabel("Time $[d]$")
ax1.set_ylabel("Transpiration rate $[kg d^{-1}]$")

ax2 = ax1.twinx()
ax2.plot(d[:, 0] / (24 * 3600), d[:, 4], 'r--')
ax2.plot(d2[:, 0] / (24 * 3600), d2[:, 4], 'g--')
ax2.plot(d3[:, 0] / (24 * 3600), d3[:, 4], 'b--')
ax2.legend(['pressure P1', 'pressure P2', 'pressure P3'], loc = 'upper right')
ax2.set_ylabel("Pressure at root collar (Pa)")

i = np.argmax(d[:, 1] == d[:, 3])
print("stress after ", d[i, 0] / (24 * 3600), "days")
i = np.argmax(d2[:, 1] == d2[:, 3])
print("stress after ", d[i, 0] / (24 * 3600), "days")
i = np.argmax(d3[:, 1] == d3[:, 3])
print("stress after ", d[i, 0] / (24 * 3600), "days")

mid1 = np.argmax(d[:, 0] / (24 * 3600) > 0.5)
mid2 = np.argmax(d2[:, 0] / (24 * 3600) > 0.5)
mid3 = np.argmax(d3[:, 0] / (24 * 3600) > 0.5)
print("transpiration during stress", d[mid1, 1], "kg/s at", d[mid1, 4], "Pa, p-crit ", d[mid1, 5], "Pa, max trans", d[mid1, 3])
print("transpiration during stress", d2[mid2, 1], "kg/s at", d2[mid2, 4], "Pa, p-crit ", d2[mid2, 5], "Pa, max trans", d2[mid2, 3])
print("transpiration during stress", d3[mid3, 1], "kg/s at", d3[mid3, 4], "Pa, p-crit ", d3[mid3, 5], "Pa, max trans", d3[mid3, 3])

plt.show()

# for 1 hour time step (2 dist)
# stress after  0.384454861111 days
# stress after  0.134454861111 days
# stress after  0.134454861111 days
# transpiration during stress 1.16697e-05 kg/s at -881000.0 Pa, p-crit  0.468517 Pa, max trans 1.16697e-05
# transpiration during stress 1.16981e-06 kg/s at -881000.0 Pa, p-crit  0.0752097 Pa, max trans 1.16981e-06
# transpiration during stress 1.07872e-06 kg/s at -881000.0 Pa, p-crit  0.0692282 Pa, max trans 1.07872e-06

# for 1 hour time step (1 dist)
# stress after  0.384454861111 days
# stress after  0.134454861111 days
# stress after  0.134454861111 days
# transpiration during stress 1.16697e-05 kg/s at -881000.0 Pa, p-crit  0.234259 Pa, max trans 1.16697e-05
# transpiration during stress 1.16981e-06 kg/s at -881000.0 Pa, p-crit  0.0376049 Pa, max trans 1.16981e-06
# transpiration during stress 1.07872e-06 kg/s at -881000.0 Pa, p-crit  0.0346141 Pa, max trans 1.07872e-06

# for 1/10 hour time step (2 dist)
# stress after  0.370833333333 days
# stress after  0.0 days
# stress after  0.0916666666667 days
# transpiration during stress 1.16697e-05 kg/s at -881000.0 Pa, p-crit  0.468517 Pa, max trans 1.16697e-05
# transpiration during stress 0.0 kg/s at -10132.0 Pa, p-crit  870868.0 Pa, max trans 13.5454
# transpiration during stress 1.07872e-06 kg/s at -881000.0 Pa, p-crit  0.0692282 Pa, max trans 1.07872e-06

# days = 1
# t_ = np.linspace(0, days, 6 * 24 * days)
# y_ = np.sin(t_ * 2.*pi - 0.5 * pi) * trans + trans
# ax1.plot(t_, y_ * (24 * 3600), 'k')
