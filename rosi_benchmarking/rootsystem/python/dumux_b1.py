#
# Compares the dumux solution of 1d single root to its analytical solution
#
# solves Benchmark 1,
# additionally computes case without gravitation, and with transpiration BC
#
# D. Leitner, 2019
#

import os
import matplotlib.pyplot as plt
from vtk_tools import *
from math import *
import van_genuchten as vg

g = 9.81  # gravitational acceleration (m/s^2)
rho = 1.e3  # density of water, (kg/m^3)
ref = 1.e5


def toPa(ph):  # cm pressure head to Pascal (kg/ (m s^2))
    return ref + ph / 100. * rho * g


def toHead(pa):  # Pascal (kg/ (m s^2)) to cm pressure head
    return (pa - ref) * 100 / rho / g


# Parameters
L = 0.5  # length of single straight root (m)
a = 2.e-3  # radius (m)
kr = 2.e-9  # radial conductivity per root type (m^2 s / kg)
kz = 5.e-13  # axial conductivity (m^5 s / kg) (mal rho ergibt die alten einheiten)
p0 = toPa(-1000)  # dircichlet bc at top (ćm)
p_s = toPa(-200)  # static soil pressure (cm)
t0 = -2e-4  # kg /s =) 17.28 kg /d
c = 2 * a * pi * kr / kz

# Boundary conditions
AA = np.array([[1, 1], [sqrt(c) * exp(sqrt(c) * (-L)), -sqrt(c) * exp(-sqrt(c) * (-L))] ])  # dirichlet top, neumann bot
bb = np.array([p0 - p_s, -rho * g])  # benchmark 1
bb2 = np.array([p0 - p_s, 0])  # neglecting gravitation

AA3 = np.array([[sqrt(c), -sqrt(c)], [sqrt(c) * exp(sqrt(c) * (-L)), -sqrt(c) * exp(-sqrt(c) * (-L))] ])  # neumann top, neumann bot
bb3 = np.array([-rho * g + t0 / 1000 / kz, -rho * g])  # transpiration as BC

d = np.linalg.solve(AA, bb)  # compute constants d_1 and d_2 from bc
d2 = np.linalg.solve(AA, bb2)  # compute constants d_1 and d_2 from bc
d3 = np.linalg.solve(AA3, bb3)  # compute constants d_1 and d_2 from bc

# Analytical solution
p_r = lambda z: toHead(p_s + d[0] * exp(sqrt(c) * z) + d[1] * exp(-sqrt(c) * z))
p_r2 = lambda z: toHead(p_s + d2[0] * exp(sqrt(c) * z) + d2[1] * exp(-sqrt(c) * z))  # neglecting gravitation
p_r3 = lambda z: toHead(p_s + d3[0] * exp(sqrt(c) * z) + d3[1] * exp(-sqrt(c) * z))  # neglecting gravitation

# Prepare plot
za_ = np.linspace(0, -L, 100)
pr = list(map(p_r, za_))
pr2 = list(map(p_r2, za_))
pr3 = list(map(p_r3, za_))

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/rosi_benchmarking/rootsystem")

# run dumux
os.system("./rootsystem input/b1.input")
p_ = read1D_vtp_data("benchmark1-00001.vtp", False)  # !!!! Box = False, CCTpfa = True
os.system("./rootsystem input/b1.input -RootSystem.Grid.File grids/singlerootH.dgf -Problem.Name benchmark1b")
p2_ = read1D_vtp_data("benchmark1b-00001.vtp", False)  # !!!! Box = False, CCTpfa = True
os.system("./rootsystem input/b1.input -RootSystem.Collar.Transpiration 17.28 -Problem.Name benchmark1c")
p3_ = read1D_vtp_data("benchmark1c-00001.vtp", False)  # !!!! Box = False, CCTpfa = True

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

# Plausibility checks
at1 = np.loadtxt("benchmark1_actual_transpiration.txt", delimiter = ',')
at2 = np.loadtxt("benchmark1b_actual_transpiration.txt", delimiter = ',')
at3 = np.loadtxt("benchmark1c_actual_transpiration.txt", delimiter = ',')

# Format of txt file:
# time_, lastActualTrans_, lastTrans_, lastMaxTrans_, p, dp, sol[0], sol[1], trans
# 0    , 1               , 2         , 3            , 4, 5,  6,      7,      8
#

print("Analytic transpiration rate (kg/s):")
print("=", 1000 * kz * ((d[0] - d[1]) * sqrt(c) + rho * g))
print("=", 1000 * kz * ((d2[0] - d2[1]) * sqrt(c) + rho * g))
print("=", 1000 * kz * ((d3[0] - d3[1]) * sqrt(c) + rho * g))
print("Numeric transpiration rate (kg/s):")
# print(at1[2])
# print(at2[2])
# print(at3[2])

print("\nAnalytic pressure at top (cm)")
print(p_r(0))
print(p_r2(0))
print(p_r3(0))
print("Numeric pressure (cm)")
# print(toHead(at1[6]))
# print(toHead(at2[6]))
# print(toHead(at3[6]))

# save benchmark 1
np.savetxt("dumux_b1", np.vstack((z_, h_)), delimiter = ',')

if __name__ == "__main__":
    plt.show()
