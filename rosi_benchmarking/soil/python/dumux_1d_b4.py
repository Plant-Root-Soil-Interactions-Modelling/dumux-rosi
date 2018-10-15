# Benchmark 4 (in 1D)
#
# compares the dumux solution to the analytical solution (Figure 5abcd Vanderborght et al 2005)
#
# D. Leitner, 2018
#

import os
import matplotlib.pyplot as plt
from analytic_b4 import *
from vtk_tools import *
from van_genuchten import *
import numpy as np
from math import *
from scipy.interpolate import interp1d

sand = Parameters(0.045, 0.43, 0.15, 3, 1.1574e-04 * 100 * 3600 * 24)
loam = Parameters(0.08, 0.43, 0.04, 1.6, 5.7870e-06 * 100 * 3600 * 24)
clay = Parameters(0.1, 0.4, 0.01, 1.1, 1.1574e-06 * 100 * 3600 * 24)

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/rosi_benchmarking/soil")

# # run dumux
# os.system("./richards1d benchmarks_1d/b4a.input")
# os.system("./richards1d benchmarks_1d/b4b.input")
# os.system("./richards1d benchmarks_1d/b4c.input")
# os.system("./richards1d benchmarks_1d/b4d.input")
# os.system("./richards1d benchmarks_1d/b4a.input -Grid.Cells 400 -Problem.Name benchmark1d_4a_hr")
# os.system("./richards1d benchmarks_1d/b4b.input -Grid.Cells 400 -Problem.Name benchmark1d_4b_hr")
# os.system("./richards1d benchmarks_1d/b4c.input -Grid.Cells 400 -Problem.Name benchmark1d_4c_hr")
# os.system("./richards1d benchmarks_1d/b4d.input -Grid.Cells 400 -Problem.Name benchmark1d_4d_hr")

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
lt = ["r:", "r:", "r:", "r:", "r-", "r-", "r-", "r-"]
for i in range(0, len(t)):
        axis[i].plot(t[i], y[i], lt[i])

np.savetxt("dumux1d", np.vstack(ex), delimiter = ",")

plt.show()

