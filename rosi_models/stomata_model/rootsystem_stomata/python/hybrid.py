#
# Felicien hybrid story
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

# print("critical pressure in cm: ", toHead(-1.4e6))


# Parameters
L = 0.5  # length of single straight root (m)
a = 2.e-3  # radius (m)
kz = 0.173  # axial conductivity [cm^4/hPa/day] similar (cm^3 / day)
kz = 1e-6 * kz / (rho * g) / (24 * 3600)
kr = 2.6e-3  # radial conductivity [cm/hPa/day] similar (1 / day)
kr = kr / (rho * g) / (24 * 3600)

p1 = toPa(-1000)  # dircichlet bc at top (ćm)
p2 = toPa(-500)  # dircichlet bc at bot (ćm)
p_s = toPa(-200)  # static soil pressure (cm)
c = 2 * a * pi * kr / kz

# Boundary conditions
AA = np.array([[1, 1], [exp(sqrt(c) * L), exp(-sqrt(c) * L)] ])  # dirichlet top, neumann bot
bb = np.array([p1 - p_s, p2 - p_s])  # benchmark 1
d = np.linalg.solve(AA, bb)  # compute constants d_1 and d_2 from bc

# Analytical solution
p_r = lambda z: toHead(p_s + d[0] * exp(sqrt(c) * z) + d[1] * exp(-sqrt(c) * z))
za_ = np.linspace(0, -L, 100)
pr = list(map(p_r, za_))
plt.plot(pr, za_, "b")
# plt.show()

print("Matrix")
print(AA)

print("\nDet")
d = exp(-sqrt(c) * L) - exp(sqrt(c) * L)
print(np.linalg.det(AA), d)

print("\nInv")
AAi = np.array([[1, 1], [exp(sqrt(c) * L), exp(-sqrt(c) * L)] ])
print(np.linalg.inv(AA))
print(AAi)
