#
# Wet top scenario for 3 phenotypes in static soil
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
# t0 = time.time()
# threads_ = []
# threads_.append(myThread("./rootsystem input/swtop.input -RootSystem.Grid.File grids/RootSys/Reference/RootSys1.dgf"))
# threads_.append(myThread("./rootsystem input/swtop.input -RootSystem.Grid.File grids/RootSys/Genotype_laterals/RootSys1.dgf -Problem.Name swtop_b"))
# threads_.append(myThread("./rootsystem input/swtop.input -RootSystem.Grid.File grids/RootSys/Genotype_volume/RootSys1.dgf -Problem.Name swtop_c"))
# for t in threads_:  # start threads
#     t.start()
# for t in threads_:  # and for all of them to finish
#      t.join()
# print("elapsed time is ", time.time() - t0)

with open("swtop_actual_transpiration.txt", 'r') as f:
    d = np.loadtxt(f, delimiter = ',')
with open("swtop_b_actual_transpiration.txt", 'r') as f:
    d2 = np.loadtxt(f, delimiter = ',')
with open("swtop_c_actual_transpiration.txt", 'r') as f:
    d3 = np.loadtxt(f, delimiter = ',')

# Format of txt file:
# time_, lastActualTrans_, lastTrans_, lastMaxTrans_, p, dp, sol[0], sol[1], trans
# 0    , 1               , 2         , 3            , 4, 5,  6,      7,      8

# Plot collar transpiration & pressure
fig, ax1 = plt.subplots()

ax1.plot(d[:, 0] / (24 * 3600), d[:, 2] * (24 * 3600) / (.75 * .15), 'k')  # potential transpiration
ax1.plot(d[:, 0] / (24 * 3600), d[:, 8] * (24 * 3600) / (.75 * .15) * 1000, 'r')  # reference, actual transpiration
ax1.plot(d[:, 0] / (24 * 3600), d[:, 3] * (24 * 3600) / (.75 * .15), 'r:')  # reference, maximal transpiration
ax1.plot(d2[:, 0] / (24 * 3600), d2[:, 8] * (24 * 3600) / (.75 * .15) * 1000, 'g')  # lateral
ax1.plot(d2[:, 0] / (24 * 3600), d2[:, 3] * (24 * 3600) / (.75 * .15), 'g:')  # lateral
ax1.plot(d3[:, 0] / (24 * 3600), d3[:, 8] * (24 * 3600) / (.75 * .15) * 1000, 'b')  # volume
ax1.plot(d3[:, 0] / (24 * 3600), d3[:, 3] * (24 * 3600) / (.75 * .15), 'b:')  # volume
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
print("transpiration during stress", d[mid1, 8], "kg/s at", d[mid1, 4], "Pa, p-crit ", d[mid1, 5], "Pa, max trans", d[mid1, 3])
print("transpiration during stress", d2[mid2, 8], "kg/s at", d2[mid2, 4], "Pa, p-crit ", d2[mid2, 5], "Pa, max trans", d2[mid2, 3])
print("transpiration during stress", d3[mid3, 8], "kg/s at", d3[mid3, 4], "Pa, p-crit ", d3[mid3, 5], "Pa, max trans", d3[mid3, 3])

plt.show()

# t_ = np.linspace(0, days, 6 * 24 * days)
# y_ = np.sin(t_ * 2.*pi - 0.5 * pi) * trans + trans
# ax1.plot(t_, y_ * (24 * 3600), 'k')

# Time step: 1/10 h (1 dist)
# stress after  0.358333333333 days
# stress after  0.3375 days
# stress after  0.3125 days
# transpiration during stress 1.12141295795e-05 kg/s at -880999.774887 Pa, p-crit  0.225113137509 Pa, max trans 1.12141295795e-05
# transpiration during stress 1.04586327942e-05 kg/s at -880999.772032 Pa, p-crit  0.227968272171 Pa, max trans 1.04586327942e-05
# transpiration during stress 9.52087489404e-06 kg/s at -880999.791658 Pa, p-crit  0.208342417609 Pa, max trans 9.52087489404e-06

# Time step: 1/10 h (2 dist)
# stress after  0.358333333333 days
# stress after  0.3375 days
# stress after  0.3125 days
# transpiration during stress 1.1214126155e-05 kg/s at -880999.549774 Pa, p-crit  0.450226137531 Pa, max trans 1.1214126155e-05
# transpiration during stress 1.04586296324e-05 kg/s at -880999.544064 Pa, p-crit  0.455936406506 Pa, max trans 1.04586296324e-05
# transpiration during stress 9.52087211702e-06 kg/s at -880999.583315 Pa, p-crit  0.41668471368 Pa, max trans 9.52087211702e-06

# INCLUDING SWITCH
# transpiration during stress 1.12141603037e-08 kg/s at 0.0 Pa, p-crit  881000.0 Pa, max trans 0.0
# transpiration during stress 1.04586075392e-08 kg/s at 0.0 Pa, p-crit  881000.0 Pa, max trans 0.0
# transpiration during stress 9.5208533841e-09 kg/s at 0.0 Pa, p-crit  881000.0 Pa, max trans 0.0

