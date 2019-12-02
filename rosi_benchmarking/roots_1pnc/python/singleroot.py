"""
 Single root test file for the root stomata model 
 
 We predescribe constant transpiration -2e-8 under constant soil pressure
"""

import os
import matplotlib.pyplot as plt
from vtk_tools import *
from math import *
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
plt.title("Predescribed transpiration")

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
# os.chdir("../../../build-cmake/rosi_benchmarking/roots_1pnc")
os.chdir("../../../build-cmake/rosi_benchmarking/roots_1pnc")

# run dumux
# os.system("./rootsystem_1pnc input/singleroot_stomata.input")
os.system("./rootsystem_1pnc input/singleroot_stomata.input")

# p_ = read1D_vtp_data("singleroot-00001.vtp")
p_ = read1D_vtp_data("singleroot-00001.vtp")

# benchmark predescribed transpiration
z_ = np.linspace(0, -0.5, len(p_))
h_ = vg.pa2head(p_)
plt.plot(h_, z_, "r+")
plt.plot(pr3, za_, "b")
plt.ylabel("Depth (m)")
plt.xlabel("Xylem pressure (cm)")
plt.title("Predescribed transpiration")

# save benchmark M31
z_ = np.linspace(0, -0.5, len(p_))
h_ = vg.pa2head(p_)
np.savetxt("dumux_m31", np.vstack((100 * z_, h_)), delimiter = ',')

with open("singleroot_actual_transpiration.txt", 'r') as f:
    d = np.loadtxt(f, delimiter = ',')

print()
c = 24 * 3600  #  [kg/s] -> [kg/per day]
print("potential", d[-1, 2] * c)
print("actual", d[-1, 1] * c)
print("actual", d[-1, 5] / 1000)
print("pressure", toHead(d[-1, 4]), pr3[0])  # root collar pressures do not perfectly agree

if __name__ == "__main__":
    plt.show()
