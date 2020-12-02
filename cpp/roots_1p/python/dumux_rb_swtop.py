#
# Wet top scenario for 3 phenotypes in static soil
#
import sys; sys.path.append("../../../../CPlantBox/src/python_modules/")

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


# first some conversions
def toHead(pa):  # Pascal (kg/ (m s^2)) to cm pressure head
    return (pa - ref) * 100 / rho / g


def toPa(ph):
    return ref + ph / 100. * rho * g;


print("wet", toPa(-100))
print("dry", toPa(-1000))

print("wet", toHead(-10000))
print("dry", toHead(-300000))
print("a week ", 7 * 24 * 3600)
trans = 5.33 * .75 * .15 / 86400  #  kg/s
maxtrans = 2 * trans
print("daily rate ", 5.33, "mm/day = ", trans, " kg/s, maximum ", maxtrans)  #
print("Critical collar pressure = ", toPa(-1.5e4))
print("kr0", np.array([1.8e-4, 1.8e-4, 0.6e-4, 0.6e-4, 0.18e-4, 0.18e-4 ]) * 1.e-4 / 86400)
print("kr1", np.array([1.8e-4, 1.8e-4, 0.18e-4, 0.18e-4 ]) * 1.e-4 / 86400)
print("kx0", np.array([0.01, 0.3, 0.3, 4.3, 4.3]) * 1.e-4 / 86400)
print("kx1", np.array([0.01e-3, 0.01e-3, 0.1e-3, 0.6e-3, 0.6e-3, 1.7e-3, 1.7e-3]) * 1.e-4 / 86400)

# Go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/cpp/roots_1p")

p1 = "-RootSystem.Grid.InitialT 63.5 -RootSystem.Grid.File maize_p1_zero_std"  # 63.5
p2 = "-RootSystem.Grid.InitialT 55.5 -RootSystem.Grid.File maize_p2_zero_std"
p3 = "-RootSystem.Grid.InitialT 58.5 -RootSystem.Grid.File maize_p3_zero_std"

# Run dumux
t0 = time.time()
threads_ = []
threads_.append(myThread("./rootsystem_rb input/rb_swtop.input " + p1 + " -Problem.Name rb_swtop_a"))
threads_.append(myThread("./rootsystem_rb input/rb_swtop.input " + p2 + " -Problem.Name rb_swtop_b"))
threads_.append(myThread("./rootsystem_rb input/rb_swtop.input " + p3 + " -Problem.Name rb_swtop_c"))
for t in threads_:  # start threads
    t.start()
for t in threads_:  # and for all of them to finish
     t.join()
print("elapsed time is ", time.time() - t0)

with open("rb_swtop_a_actual_transpiration.txt", 'r') as f:
    d = np.loadtxt(f, delimiter = ',')
with open("rb_swtop_b_actual_transpiration.txt", 'r') as f:
    d2 = np.loadtxt(f, delimiter = ',')
with open("rb_swtop_c_actual_transpiration.txt", 'r') as f:
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
# ax2.plot(t, d[:, 4], 'r--')
# ax2.plot(t, d2[:, 4], 'g--')
# ax2.plot(t, d3[:, 4], 'b--')
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

