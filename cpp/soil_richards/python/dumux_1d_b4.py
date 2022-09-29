# Benchmark 4 (in 1D)
#
# compares the dumux solution to the analytical solution (Figure 5abcd Vanderborght et al 2005)
#
# D. Leitner, 2018
#
import sys; sys.path.append("../../../../CPlantBox/src/python_modules/")
import sys; sys.path.append("../../../python/soil/")  # for the analytical solutions

import os
import matplotlib.pyplot as plt
from analytic_b4 import *
from vtk_tools import *
from van_genuchten import *
import numpy as np
from math import *
from scipy.interpolate import interp1d

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/cpp/soil_richards")

# run dumux
os.system("./richards1d input/b4a_1d.input")
os.system("./richards1d input/b4b_1d.input")
os.system("./richards1d input/b4c_1d.input")
os.system("./richards1d input/b4d_1d.input")
os.system("./richards1d input/b4a_1d.input -Soil.Grid.Cells 1399 -Problem.Name benchmark1d_4a_hr")
os.system("./richards1d input/b4b_1d.input -Soil.Grid.Cells 1399 -Problem.Name benchmark1d_4b_hr")
os.system("./richards1d input/b4c_1d.input -Soil.Grid.Cells 1399 -Problem.Name benchmark1d_4c_hr")
os.system("./richards1d input/b4d_1d.input -Soil.Grid.Cells 1399 -Problem.Name benchmark1d_4d_hr")

# open results
num = ['a', 'c', 'b', 'd', 'a_hr', 'c_hr', 'b_hr', 'd_hr']
t = []
y = []
ex = []
minN = inf

for i, n in enumerate(num):
    try:
        with open("benchmark1d_4" + n + ".csv", 'r') as f:
            d = np.loadtxt(f, delimiter = ',')
            t_ = d[:, 0] / (24 * 3600)
            f_ = d[:, 1] * (24 * 3600) / 1000 * 100;  # from kg/(s m²) to cm/day
            t_, I = np.unique(t_, return_index = True)
            f_ = f_[I]
            t.append(t_)
            y.append(f_)
            if (len(t_) < minN):
                minN = len(t_)
    except:
        pass

# resample data to minN points
print("Resampling to", minN, "data points")
for i in range(0, len(t)):

    t_ = t[i]
    y_ = y[i]

    inter = interp1d(t_, y_, kind = 'linear', fill_value = 'extrapolate')
    new_t = np.linspace(np.min(t_), np.max(t_), minN)
    new_y = inter(new_t)
    ex.append(new_t)
    ex.append(new_y)

# prepare plot
axis = [ax1, ax2, ax3, ax4, ax1, ax2, ax3, ax4]
lt = ["r", "r", "r", "r", "g--", "g--", "g--", "g--"]
for i in range(0, len(t)):
        axis[i].plot(t[i], y[i], lt[i])

np.savetxt("dumux1d_b4", np.vstack(ex), delimiter = ",")

ax1.set_ylim(0, 0.11)
ax2.set_ylim(0, 0.11)
ax3.set_ylim(0, 0.31)
ax4.set_ylim(0, 0.31)

# ax2.set_xlim(0, 10)
# ax2.set_xlim(0, 10)
# ax2.set_xlim(0, 10)

# plt.show()

