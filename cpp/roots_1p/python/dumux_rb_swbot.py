#
# Wet bot scenario for 3 phenotypes in static soil
#
import sys; sys.path.append("../../../../CPlantBox/src/python_modules/")

import os
import matplotlib.pyplot as plt
from vtk_tools import *
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


# Go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/cpp/roots_1p")

p1 = "-RootSystem.Grid.InitialT 63.5 -RootSystem.Grid.File maize_p1_zero_std"  # 63.5 File maize_p1_zero_std    Anagallis_femina_Leitner_2010
p2 = "-RootSystem.Grid.InitialT 55.5 -RootSystem.Grid.File maize_p2_zero_std"
p3 = "-RootSystem.Grid.InitialT 58.5 -RootSystem.Grid.File maize_p3_zero_std"

# Run dumux
t0 = time.time()
threads_ = []
threads_.append(myThread("./rootsystem_rb input/rb_swbot.input " + p1 + " -Problem.Name rb_swbot_a"))
threads_.append(myThread("./rootsystem_rb input/rb_swbot.input " + p2 + " -Problem.Name rb_swbot_b"))
threads_.append(myThread("./rootsystem_rb input/rb_swbot.input " + p3 + " -Problem.Name rb_swbot_c"))
for t in threads_:  # start threads
    t.start()
for t in threads_:  # and for all of them to finish
     t.join()
print("elapsed time is ", time.time() - t0)

# 0 time [s], 1 actual transpiration [kg/s], 2 potential transpiration [kg/s], 3 maximal transpiration [kg/s],
# 4 collar pressure [Pa], 5 calculated actual transpiration
with open("rb_swbot_a_actual_transpiration.txt", 'r') as f:
    d = np.loadtxt(f, delimiter = ',')
with open("rb_swbot_b_actual_transpiration.txt", 'r') as f:
    d2 = np.loadtxt(f, delimiter = ',')
with open("rb_swbot_c_actual_transpiration.txt", 'r') as f:
    d3 = np.loadtxt(f, delimiter = ',')

# Plot collar transpiration & pressure
fig, ax1 = plt.subplots()

c = 24 * 3600  # [kg/s] -> [kg/day]
t = d[:, 0] / (24 * 3600)
t2 = d2[:, 0] / (24 * 3600)
t3 = d3[:, 0] / (24 * 3600)

ax1.plot(t, d[:, 2] * c, 'k')  # potential transpiration  * c

ax1.plot(t, d[:, 1] * c, 'r-,')  # reference, actual transpiration
ax1.plot(t2, d2[:, 1] * c, 'g-,')  # lateral
ax1.plot(t3, d3[:, 1] * c, 'b-,')  # volume
# ax1.plot(t, d[:, 5]  , 'r:,')  # calculated actual transpiration
# print(d[:, 5] / 100000)

print(d[:, 5])

ax1.legend(['Pot trans', 'actual trans P1', 'actual trans P2', 'actual trans P3'], loc = 'upper left')
ax1.axis((0, t[-1], 0, 1.2))
ax1.set_xlabel("Time $[d]$")
ax1.set_ylabel("Transpiration rate $[kg \ d^{-1}]$")

# ax2 = ax1.twinx()
# ax2.plot(t, toHead(d[:, 4]), 'r--')
# ax2.plot(t2, toHead(d2[:, 4]), 'g--')
# ax2.plot(t3, toHead(d3[:, 4]), 'b--')
# ax2.legend(['pressure P1', 'pressure P2', 'pressure P3'], loc = 'upper right')
# ax2.set_ylabel("Pressure at root collar (cm)")

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

