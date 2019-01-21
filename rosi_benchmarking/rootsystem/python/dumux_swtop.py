#
# Wet top scenario
#

import os
import matplotlib.pyplot as plt
from vtk_tools import *
from math import *
import van_genuchten as vg

g = 9.81  # gravitational acceleration (m/s^2)
rho = 1.e3  # density of water, (kg/m^3)
ref = 1.e5


def toHead(pa):  # Pascal (kg/ (m s^2)) to cm pressure head
    return (pa - ref) * 100 / rho / g

# run dumux


# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/rosi_benchmarking/rootsystem")

print("wet", toHead(-10000))
print("dry", toHead(-300000))
print("a week ", 7 * 24 * 3600)

trans = 5.33 * .75 * .15 / 86400
mtrans = 2 * trans
print("daily rate ", 5.33, "mm/day = ", trans, " kg/s, maximum ", mtrans)  #

# run dumux
# os.system("./rootsystem input/swtop.input")  # mpirun -n 8
p_, z_ = read3D_vtp_data("swtop-00000.vtp", False)
h_ = vg.pa2head(p_)
plt.plot(p_, z_[:, 2], "r+")  # cell data
plt.ylabel("Depth (m)")
plt.xlabel("Xylem pressure (cm)")
plt.show()

# with open("swtop_actual_transpiration.txt", 'r') as f:
#     d = np.loadtxt(f, delimiter = ',')
#
# t_ = np.linspace(0, 7, 6 * 24 * days)
# y_ = np.sin(t_ * 2.*pi - 0.5 * pi) * trans + trans
#
# plt.plot(t_, y_ * (24 * 3600), 'k')
# plt.plot(d[:, 0] / (24 * 3600), d[:, 1] * (24 * 3600), 'r')
# plt.xlabel("Time $[d]$")
# plt.ylabel("Transpiration rate $[mm d^{-1}]$")
# plt.show()

