# Run Benchmakrk C11

import os
import matplotlib.pyplot as plt
from vtk_tools import *
import van_genuchten as vg
import math

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/rosi_benchmarking/coupled")

soiltype = 2  # sand, loam, clay
q_root = 0.1
soils = ["Sand", "Loam", "Clay"]
trans = q_root * (2 * 0.02 * math.pi) * 1 * 1.e-6 * 1000.  # cm/day -> kg/day
print(trans)
os.system("./coupled input/benchmarkC11.input -Soil.Layer.Number {} -RootSystem.Collar.Transpiration {}"
          .format(soiltype, trans))

# Figure
leg_str = [];
for i in range(0, 7):
    p_, y_ = read3D_vtp_data_line("benchmarkC11-0000{}.vtu".format(i), False)
    h1_ = vg.pa2head(p_)
    plt.plot(np.array(y_) * 100, h1_)
    leg_str.append("{}".format(i * 10) + " days")

plt.legend(leg_str)
plt.ylabel('$\psi$ (cm)')
plt.xlabel('y axis (cm)')
plt.title(soils[soiltype] + ", q_{root}=" + str(q_root) + "d")
plt.show()

