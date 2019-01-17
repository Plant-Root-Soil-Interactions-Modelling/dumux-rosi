#
# compares the dumux solution of 1d root model to its analytical solution
#
# D. Leitner, 2018
#

import os
import matplotlib.pyplot as plt
from vtk_tools import *
from math import *
import van_genuchten as vg

#
# Analytical solution
#

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

# Analytical solution
c = 2 * a * pi * kr / kz
p_r = lambda z: toHead(p_s + d[0] * exp(sqrt(c) * z) + d[1] * exp(-sqrt(c) * z))

# Boundary conditions
AA = np.array([[1, 1], [sqrt(c) * exp(sqrt(c) * (-L)), -sqrt(c) * exp(-sqrt(c) * (-L))] ])  # dirichlet top, neumann bot
bb = np.array([p0 - p_s, -rho * g])  #
d = np.linalg.solve(AA, bb)  # compute constants d_1 and d_2 from bc

print("analytic transpiration rate", 1000 * kz * (d[0] - d[1]) * sqrt(c), "kg/s")
print("= ", 1000 * kz * (d[0] - d[1]) * sqrt(c) / .75 / .15 * 86400, "mm day-1")

# Prepare plot
za_ = np.linspace(0, -L, 100)
pr = list(map(p_r, za_))

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/rosi_benchmarking/rootsystem")

# run dumux
os.system("./rootsystem input/b1.input")
p_ = read1D_vtp_data("benchmark1-00001.vtp", False)  # !!!! Box = False, CCTpfa = True
z_ = np.linspace(0, -0.5, len(p_))
h_ = vg.pa2head(p_)
plt.plot(h_, z_, "r+")
plt.plot(pr, za_)
plt.ylabel("Depth (m)")
plt.xlabel("Xylem pressure (cm)")
np.savetxt("dumux_b1", np.vstack((z_, h_)), delimiter = ',')
plt.show()

