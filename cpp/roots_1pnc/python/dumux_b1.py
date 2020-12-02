"""
TODO not working

 Benchmark M31

 Compares the dumux solution of 1d single root to its analytical solution

 solves Benchmark M31,
 additionally computes case without gravitation, and with transpiration BC

 D. Leitner, 2019
"""
import sys; sys.path.append("../../../../CPlantBox/src/python_modules/")

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

p0 = toPa(-1000)  # dircichlet bc at top (ćm)
p_s = toPa(-200)  # static soil pressure (cm)
t0 = -2e-8  # kg / s
trans = -t0 * 24 * 3600  # kg /day
print("tranpsiration ", trans, "[kg/day]")

# Boundary conditions
c = 2 * a * pi * kr / kz
AA = np.array([[1, 1], [sqrt(c) * exp(sqrt(c) * (-L)), -sqrt(c) * exp(-sqrt(c) * (-L))] ])  # dirichlet top, neumann bot
bb = np.array([p0 - p_s, -rho * g])  # benchmark 1
bb2 = np.array([p0 - p_s, 0])  # neglecting gravitation

AA3 = np.array([[sqrt(c), -sqrt(c)], [sqrt(c) * exp(sqrt(c) * (-L)), -sqrt(c) * exp(-sqrt(c) * (-L))] ])  # neumann top, neumann bot
bb3 = np.array([-rho * g + t0 / rho / kz, -rho * g])  # transpiration as BC

# Integration constants
d = np.linalg.solve(AA, bb)  # compute constants d_1 and d_2 from bc
d2 = np.linalg.solve(AA, bb2)  # compute constants d_1 and d_2 from bc
d3 = np.linalg.solve(AA3, bb3)  # compute constants d_1 and d_2 from bc

# Analytical solution
p_r = lambda z: toHead(p_s + d[0] * exp(sqrt(c) * z) + d[1] * exp(-sqrt(c) * z))
p_r2 = lambda z: toHead(p_s + d2[0] * exp(sqrt(c) * z) + d2[1] * exp(-sqrt(c) * z))  # neglecting gravitation
p_r3 = lambda z: toHead(p_s + d3[0] * exp(sqrt(c) * z) + d3[1] * exp(-sqrt(c) * z))  # transpiration

# Prepare plot
za_ = np.linspace(0, -L, 100)
pr = list(map(p_r, za_))
pr2 = list(map(p_r2, za_))
pr3 = list(map(p_r3, za_))

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/cpp/roots_1pnc")

# delete old results
os.system("rm benchmark1-00001.vtp")
# os.system("rm benchmark1b-00001.vtp")
# os.system("rm benchmark1c-00001.vtp")

# run dumux
os.system("./rootsystem_stomata2 input/b1.input")
p_ = read1D_vtp_data("benchmark1-00001.vtp")
# os.system("./rootsystem_stomata2 input/b1.input -RootSystem.Grid.File grids/singlerootH.dgf -Problem.Name benchmark1b")
# p2_ = read1D_vtp_data("benchmark1b-00001.vtp")
# os.system("./rootsystem_stomata2 input/b1_trans.input -RootSystem.Collar.Transpiration {} -Problem.Name benchmark1c".format(trans))
# p3_ = read1D_vtp_data("benchmark1c-00001.vtp")

# plot
fig, (ax1, ax2, ax3) = plt.subplots(1, 3)

# benchmark 1
z_ = np.linspace(0, -0.5, len(p_))
h_ = vg.pa2head(p_)
ax1.plot(h_, z_, "r+")
ax1.plot(pr, za_, "b")
ax1.set_ylabel("Depth (m)")
ax1.set_xlabel("Xylem pressure (cm)")
ax1.set_title("Benchmark 1")

# benchmark without gravity
h_ = vg.pa2head(p2_)
ax2.plot(h_, z_, "r+")
ax2.plot(pr, za_, "c")
ax2.plot(pr2, za_, "b")
ax2.set_ylabel("Depth (m)")
ax2.set_xlabel("Xylem pressure (cm)")
ax2.set_title("Benchmark 1 - neglecting gravitation")

# benchmark predescribed transpiration
h_ = vg.pa2head(p3_)
ax3.plot(h_, z_, "r+")
ax3.plot(pr3, za_, "b")
ax3.set_ylabel("Depth (m)")
ax3.set_xlabel("Xylem pressure (cm)")
ax3.set_title("Predescribed transpiration")

# save benchmark M31
z_ = np.linspace(0, -0.5, len(p_))
h_ = vg.pa2head(p_)
np.savetxt("dumux_m31", np.vstack((100 * z_, h_)), delimiter = ',')

with open("benchmark1c_actual_transpiration.txt", 'r') as f:
    d = np.loadtxt(f, delimiter = ',')

print()
c = 24 * 3600  #  [kg/s] -> [kg/per day]
print("potential", d[-1, 2] * c)
print("actual", d[-1, 1] * c)
print("actual", d[-1, 5] / 1000)
print("pressure", toHead(d[-1, 4]), pr3[0])  # root collar pressures do not perfectly agree

if __name__ == "__main__":
    plt.show()
