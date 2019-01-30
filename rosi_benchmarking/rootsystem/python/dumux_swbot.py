#
# Bot top scenario
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


# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/rosi_benchmarking/rootsystem")

print("wet", toHead(-10000))
print("dry", toHead(-300000))
print("a week ", 7 * 24 * 3600)

trans = 5.33 * .75 * .15 / 86400
maxtrans = 2 * trans
print("daily rate ", 5.33, "mm/day = ", trans, " kg/s, maximum ", maxtrans)  #

# run dumux
os.system("./rootsystem input/swbot.input -RootSystem.Collar.Transpiration \"0 1.38802083e-05\"")  # mpirun -n 8 1.38802083e-05

with open("swbot_actual_transpiration.txt", 'r') as f:
    d = np.loadtxt(f, delimiter = ',')

days = 1
t_ = np.linspace(0, days, 6 * 24 * days)
y_ = np.sin(t_ * 2.*pi - 0.5 * pi) * trans + trans

# Plot collar transpiration & pressure
fig, ax1 = plt.subplots()

ax1.plot(t_, y_ * (24 * 3600), 'k')
ax1.plot(d[:, 0] / (24 * 3600), d[:, 1] * (24 * 3600), 'b')
# ax1.plot(d[:, 0] / (24 * 3600), d[:, 2] * (24 * 3600), 'r')
ax1.set_xlabel("Time $[d]$")
ax1.set_ylabel("Transpiration rate $[kg d^{-1}]$")

ax2 = ax1.twinx()
ax2.plot(d[:, 0] / (24 * 3600), d[:, 3], 'g')
ax2.set_ylabel("Pressure at root collar (cm)")
plt.show()

# # plot
# p_, z_ = read3D_vtp_data("swbot-00000.vtp", False)
# h_ = vg.pa2head(p_)
# plt.plot(p_, z_[:, 2], "r+")  # cell data
# plt.ylabel("Depth (m)")
# plt.xlabel("Xylem pressure (cm)")
# plt.show()

