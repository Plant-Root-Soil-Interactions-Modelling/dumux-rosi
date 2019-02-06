#
# Wet bot scenario for 3 phenotypes in static soil
#
import os
import matplotlib.pyplot as plt
from vtk_tools import *
from math import *
import van_genuchten as vg
import threading
import time

g = 9.81  # gravitational acceleration (m/s^2)
rho = 1.e3  # density of water, (kg/m^3)
ref = 1.e5


class myThread(threading.Thread):

    def __init__(self, name):
        threading.Thread.__init__(self)
        self.name = name

    def run(self):
        print("Starting:" + self.name)
        os.system(self.name)


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

with open("swbot_actual_transpiration.txt", 'r') as f:
    d = np.loadtxt(f, delimiter = ',')
with open("swbot_b_actual_transpiration.txt", 'r') as f:
    d2 = np.loadtxt(f, delimiter = ',')
with open("swbot_c_actual_transpiration.txt", 'r') as f:
    d3 = np.loadtxt(f, delimiter = ',')

# Format of txt file:
# time_, lastActualTrans_, lastTrans_, lastMaxTrans_, p, dp, sol[0], sol[1], trans
# 0    , 1               , 2         , 3            , 4, 5,  6,      7,      8
#

# Plot collar transpiration & pressure
fig, ax1 = plt.subplots()

ax1.plot(d[:, 0] / (24 * 3600), d[:, 2] * (24 * 3600) / (.75 * .15), 'k')  # potential transpiration
ax1.plot(d[:, 0] / (24 * 3600), d[:, 1] * (24 * 3600) / (.75 * .15), 'r')  # reference, actual transpiration
ax1.plot(d[:, 0] / (24 * 3600), d[:, 3] * (24 * 3600) / (.75 * .15), 'r:')  # reference, maximal transpiration
ax1.plot(d2[:, 0] / (24 * 3600), d2[:, 1] * (24 * 3600) / (.75 * .15), 'g')  # lateral
ax1.plot(d2[:, 0] / (24 * 3600), d2[:, 3] * (24 * 3600) / (.75 * .15), 'g:')  # lateral
ax1.plot(d3[:, 0] / (24 * 3600), d3[:, 1] * (24 * 3600) / (.75 * .15), 'b')  # volume
ax1.plot(d3[:, 0] / (24 * 3600), d3[:, 1] * (24 * 3600) / (.75 * .15), 'b:')  # volume
ax1.legend(['Pot trans', 'actual trans P1', 'max trans P1', 'actual trans P2', 'max trans P2', 'actual trans P3', 'max trans P3'], loc = 'upper left')
ax1.axis((0, d[-1, 0] / (24 * 3600), 0, 12))
ax1.set_xlabel("Time $[d]$")
ax1.set_ylabel("Transpiration rate $[mm \ d^{-1}]$")

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

# for 1 hour time step (1 dist)
# stress after  0.345833333333 days
# stress after  0.3125 days
# stress after  0.316666666667 days
# transpiration during stress 1.07743732989e-05 kg/s at -880999.783715 Pa, p-crit  0.216285442468 Pa, max trans 1.07743732989e-05
# transpiration during stress 9.5668974017e-06 kg/s at -880999.791469 Pa, p-crit  0.20853099192 Pa, max trans 9.5668974017e-06
# transpiration during stress 9.72600943435e-06 kg/s at -880999.787169 Pa, p-crit  0.212831314537 Pa, max trans 9.72600943435e-06

# for 1 hour time step (2 dist)
# stress after  0.345833333333 days
# stress after  0.3125 days
# stress after  0.316666666667 days
# transpiration during stress 1.07743700078e-05 kg/s at -880999.567429 Pa, p-crit  0.432570752804 Pa, max trans 1.07743700078e-05
# transpiration during stress 9.5668945123e-06 kg/s at -880999.582938 Pa, p-crit  0.417061857879 Pa, max trans 9.5668945123e-06
# transpiration during stress 9.72600659614e-06 kg/s at -880999.574337 Pa, p-crit  0.425662504858 Pa, max trans 9.72600659614e-06

# days = 1
# t_ = np.linspace(0, days, 6 * 24 * days)
# y_ = np.sin(t_ * 2.*pi - 0.5 * pi) * trans + trans
# ax1.plot(t_, y_ * (24 * 3600), 'k')
