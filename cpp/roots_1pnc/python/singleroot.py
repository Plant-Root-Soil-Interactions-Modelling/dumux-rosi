"""
 Single root test file for the root stomata model 
 
 We predescribe constant transpiration -2e-8 under constant soil pressure
"""
import sys; sys.path.append("../../../../CPlantBox/src/python_modules/")

import os
from math import *
import numpy as np
import matplotlib.pyplot as plt
from vtk_tools import *
import van_genuchten as vg

g = 9.81  # gravitational acceleration (m/s^2)
rho = 1.e3  # density of water, (kg/m^3)
ref = 1.e5  # Pa


def toPa(ph):  # cm pressure head to Pascal (kg/ (m s^2))
    return ref + ph / 100. * rho * g


def toHead(pa):  # Pascal (kg/ (m s^2)) to cm pressure head
    return (pa - ref) * 100 / rho / g


L = 0.5  # length of single straight root (m)
a = 2.e-3  # radius (m)
kz = 4.32e-2  # axial conductivity [cm^4/hPa/day] similar (cm^3 / day)
kz = 1e-6 * kz / (rho * g) / (24 * 3600)
kr = 1.728e-4  # radial conductivity [cm/hPa/day] similar (1 / day)
kr = kr / (rho * g) / (24 * 3600)

p0 = toPa(-1000)  # dircichlet bc at top (ćm)
p_s = toPa(-200)  # static soil pressure (cm)
t0 = -2e-8  # kg / s
trans = -t0 * 24 * 3600  # kg /day
print(trans)

# Boundary conditions
c = 2 * a * pi * kr / kz
AA3 = np.array([[sqrt(c), -sqrt(c)], [sqrt(c) * exp(sqrt(c) * (-L)), -sqrt(c) * exp(-sqrt(c) * (-L))] ])  # neumann top, neumann bot
bb3 = np.array([-rho * g + t0 / rho / kz, -rho * g])  # transpiration as BC

# Integration constants
d3 = np.linalg.solve(AA3, bb3)  # compute constants d_1 and d_2 from bc

p_r3 = lambda z: toHead(p_s + d3[0] * exp(sqrt(c) * z) + d3[1] * exp(-sqrt(c) * z))  # transpiration

# Prepare plot
za_ = np.linspace(0, -L, 100)
pr3 = list(map(p_r3, za_))

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/cpp/roots_1pnc")

# run dumux
os.system("./rootsystem_stomata2 input/singleroot_stomata.input")

""" benchmark pressure head in single root """
p_, z_ = read3D_vtp_data("singleroot-00001.vtp")
h_ = vg.pa2head(p_)
plt.plot(h_, z_[:, 2], "r+")
plt.plot(pr3, za_, "b")
plt.ylabel("Depth (m)")
plt.xlabel("Xylem pressure (cm)")

#      * 0 time [s], 1 actual transpiration [kg/s], 2 potential transpiration [kg/s], 3 maximal transpiration [kg/s],
#      * 4 collar pressure [Pa], 5 calculated actual transpiration [cm^3/day], 6 simtime [s], 7 hormone leaf mass [kg],
#      * 8 hormone collar flow rate [kg/s], 9 hormone root system mass [kg] , 10 hormone source rate [kg/s]
with open("singleroot_actual_transpiration.txt", 'r') as f:
    d = np.loadtxt(f, delimiter = ',')
c = 24 * 3600  # s / day

# """ Plot transpiration """
# plt.plot(d[:, 0] / c, 1000 * d[:, 2] * c, 'k')  # potential transpiration
# plt.plot(d[:, 0] / c, 1000 * d[:, 1] * c, 'r-,')  # actual transpiration
# plt.xlabel("time (days)")
# plt.ylabel("transpiration (g/day)")
# plt.legend(["potential", "actual"])
# plt.title("Water transpiration")

""" Plot hormone rate and mass """
# fig, [ax1, ax2] = plt.subplots(1, 2)
#
# ax1.plot(d[:, 0] / c, 1000 * d[:, 8] * c, "r")
# ax1.plot(d[:, 0] / c, 1000 * d[:, 10] * c, "b")
# ax1.set_ylabel("mass rate (g/day)")
# ax1.set_xlabel("time (days)")
# ax1.set_title("Hormone production rate")
# ax1.legend(["leaf rate", "root system rate"])
#
# ax2.plot(d[:, 0] / c, 1000 * d[:, 7], "r")
# ax2.plot(d[:, 0] / c, 1000 * d[:, 9], "b")
# ax2.set_ylabel("mass (g)")
# ax2.set_xlabel("time (days)")
# ax2.set_title("Hormone mass")
# ax2.legend(["leaf mass", "root system mass"])

if __name__ == "__main__":
    plt.show()
