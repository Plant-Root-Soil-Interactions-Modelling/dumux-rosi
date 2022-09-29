''' Run Benchmakrk C11 (classic sink)'''

""" todo: check results, make a nicer plot """
import sys; sys.path.append("../../../../CPlantBox/src/python_modules/")

import os
import matplotlib.pyplot as plt
from vtk_tools import *
import van_genuchten as vg
import math

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/cpp/coupled_1p_richards")

soils = ["Sand", "Loam", "Clay"]


def simulate(soiltype, q_root, simtime, checktimes, maxtimestep):
    trans = q_root * ((2 * 0.2 * math.pi) * 1) * 1.e-6 * 1000.  # cm/day -> kg/day
    print("Transpiration", trans, "kg/day")
    os.system("./coupled input/benchmarkC11.input -TimeLoop.TEnd {} -TimeLoop.PeriodicCheckTimes {} -TimeLoop.MaxTimeStepSize {} ".format(simtime, checktimes, maxtimestep)
              +"-Problem.Name {} -Soil.Layer.Number {} -RootSystem.Collar.Transpiration {}".format("benchmarkC11" + soils[soiltype], soiltype + 1, trans))


def plot_number(ax1, i, tend, soiltype, q_root, export):

    print("benchmarkC11" + soils[soiltype] + "-000{0:02d}.vtu".format(i))
    p_, y_ = read3D_data("benchmarkC11" + soils[soiltype] + "-000{0:02d}".format(i), 1, 2)
    h1_ = vg.pa2head(p_) 
    ax1.plot(np.array(y_[:,0]) * 100, h1_)
    export.append(np.array(y_[:,0]) * 100)
    export.append(h1_)
    ax1.legend([str(tend) + " days"])
    ax1.set_ylabel('$water potential$ (cm)')
    ax1.set_xlabel('r (cm)')
    ax1.set_title(soils[soiltype] + ", q_{root}=" + str(q_root) + "cm/d")
    #ax1.axis([0., 0.5, -5000, 0.])

def plot_actual(ax1, soiltype, q_root):
    with open("benchmarkC11" + soils[soiltype] + "_actual_transpiration.txt", 'r') as f:
        d = np.loadtxt(f, delimiter = ',')
    c = (24 * 3600)  #  kg/s -> cm^3 / per day
    t = d[:, 0] / (24 * 3600)
    # 0 time, 1 actual transpiration, 2 potential transpiration, 3 maximal transpiration, 4 collar pressure, 5 calculated actual transpiration
    ax1.plot(t, d[:, 2] * c, 'k')  # potential transpiration
    ax1.plot(t, d[:, 1] * c, 'r-,')  # reference, actual transpiration
    ax1.set_title(soils[soiltype] + ", q_{root}=" + str(q_root) + "d")
    ax1.legend(['Pot trans', 'actual trans'])
    ax1.set_xlabel("Time $[d]$")
    ax1.set_ylabel("Transpiration rate $[kg \ d^{-1}]$")


def stress_at(soiltype, q_root):
    with open("benchmarkC11" + soils[soiltype] + "_actual_transpiration.txt", 'r') as f:
        d = np.loadtxt(f, delimiter = ',')
    c = (24 * 3600)  #  kg/s -> cm^3 / per day
    t = d[:, 0] / (24 * 3600)
    stress_t = t[np.argwhere(d[:, 1] < d[:, 2])[0]]
    print("First stress occurs at day", stress_t)
    return float(stress_t[0])

# q_root = 0.1  # cm/day
# fig, axes = plt.subplots(1, 3, sharey = True)
# for i in range(0, 3):
#
#     soiltype = i  # sand, loam, clay
#
#     # simulate for 17 days, daily check times, max time step 1h
#     simulate(soiltype, q_root, 17 * 24 * 3600, 24 * 3600, 3600)
#
#     # find when the stress occured
#     stress_t = stress_at(soiltype, q_root)  # get stress time
#
#     n = (int)(stress_t + 0.99)
#     plot_number(axes[i], n + 1, stress_t, soiltype, q_root)
#
#     # # simulate that point in time again
#     # tend = stress_t * 24 * 3600
#     # simulate(soiltype, q_root, tend, tend, 3600)
#     # plot_number(ax1, 1, stress_t, soiltype, q_root)
#
# plt.show()


export = []

q = [0.1, 0.05]  # cm/day
times = np.array([[0.1, 0.1], [10, 21.2], [8.5, 17.5]])

fig, axes = plt.subplots(2, 3, sharey = True)

for j in range(0, 2):
    for i in range(0, 3):
        tend = times[i][j] * 24 * 3600
        soiltype = i  # sand, loam, clay
        q_root = q[j]
        simulate(soiltype, q_root, tend, tend, 3600)
        plot_number(axes[j, i], 1, times[i][j], soiltype, q_root, export)

# plt.show()

# for e in export:
#     print(e.shape)

np.savetxt("dumux_c11", np.vstack(export), delimiter = ",")

# ex.append(z_)
# ex.append(theta_)
# np.savetxt("dumux1d_b3", np.vstack(ex), delimiter = ",")
# ex.append(z_)
# ex.append(theta_)
# np.savetxt("dumux1d_b3", np.vstack(ex), delimiter = ",")
