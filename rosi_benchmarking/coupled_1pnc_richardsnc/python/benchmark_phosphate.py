''' Run Benchmakrk C11 '''

import os
import matplotlib.pyplot as plt
from vtk_tools import *
import van_genuchten as vg
import math

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/rosi_benchmarking/coupled_1pnc_richardsnc")

trans = 1.e-7 * ((2 * 0.02 * math.pi) * 1) * 1.e-6 * 1000. * (24.*3600.)  # cm/s -> kg/day
print("Transpiration", trans, "kg/day")

os.system("./coupled_1p2c input/benchmark_phosphate.input")

 # -TimeLoop.TEnd {} -TimeLoop.PeriodicCheckTimes {} -TimeLoop.MaxTimeStepSize {} ".format(simtime, checktimes, maxtimestep)
              # +"-Problem.Name {} -Soil.Layer.Number {} -RootSystem.Collar.Transpiration {}".format("benchmarkC11" + soils[soiltype], soiltype + 1, trans))

#     print("benchmarkC11" + soils[soiltype] + "-000{0:02d}.vtu".format(i))
#     p_, y_ = read3D_vtp_data_line("benchmarkC11" + soils[soiltype] + "-000{0:02d}.vtu".format(i), 1.e-4)
#     h1_ = vg.pa2head(p_)  #
#     print("h", np.min(h1_), np.max(h1_))
#     print("y", np.min(y_), np.max(y_))
#
#     ax1.plot(np.array(y_) * 100, h1_)
#     export.append(np.array(y_) * 100)
#     export.append(h1_)
#
#     ax1.legend([str(tend) + " days"])
#     ax1.set_ylabel('$water potential$ (cm)')
#     ax1.set_xlabel('r (cm)')
#     ax1.set_title(soils[soiltype] + ", q_{root}=" + str(q_root) + "cm/d")
#     ax1.axis([0., 0.5, -5000, 0.])
#
#
# def plot_actual(ax1, soiltype, q_root):
#     with open("benchmarkC11" + soils[soiltype] + "_actual_transpiration.txt", 'r') as f:
#         d = np.loadtxt(f, delimiter = ',')
#     c = (24 * 3600)  #  kg/s -> cm^3 / per day
#     t = d[:, 0] / (24 * 3600)
#     # 0 time, 1 actual transpiration, 2 potential transpiration, 3 maximal transpiration, 4 collar pressure, 5 calculated actual transpiration
#     ax1.plot(t, d[:, 2] * c, 'k')  # potential transpiration
#     ax1.plot(t, d[:, 1] * c, 'r-,')  # reference, actual transpiration
#     ax1.set_title(soils[soiltype] + ", q_{root}=" + str(q_root) + "d")
#     ax1.legend(['Pot trans', 'actual trans'])
#     ax1.set_xlabel("Time $[d]$")
#     ax1.set_ylabel("Transpiration rate $[kg \ d^{-1}]$")

