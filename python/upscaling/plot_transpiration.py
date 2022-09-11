"""
transpiration plot (one column, number of rows as number of filenames)
"""
import sys; sys.path.append("../../../CPlantBox/");  sys.path.append("../../../CPlantBox/src/python_modules")

import numpy as np
import matplotlib.pyplot as plt

from xylem_flux import sinusoidal2

# name = "singleroot"
# # str_ = ["cyl", "sra"]
# # titles = ["cylindrical", "steady rate"]
# str_ = ["sra", "ups", "agg"]
# titles = ["steady rate", "upscaling", "parallel"]  # "cylindrical"
# # str_ = ["cyl", "sra", "agg", "upp"]
# # titles = ["cylindrical", "steady rate", "parallel", "upscaled"]
# trans = 0.6 * 4  # potential transpiration cm3/day
# fnames = np.array(["transpiration_" + name + "_" + s for s in str_ ])

name = "rootsystem1d"
str_ = ["sra", "ups", "agg"]
titles = ["steady rate", "upscaled", "parallel"]
trans = 0.6 * (12 * 3)  # potential transpiration cm3/day
fnames = np.array(["transpiration_" + name + "_" + s for s in str_ ])

# names = ["rootsystem_nobasals1d_sra", "rootsystem_nobasals1d_agg"]
# str_ = ["sra", "agg"]
# titles = names
# trans = 0.6 * (12 * 3)  # potential transpiration cm3/day
# fnames = np.array(["transpiration_" + name for name in names])

path = "results/"

SMALL_SIZE = 16
MEDIUM_SIZE = 16
BIGGER_SIZE = 16
plt.rc('font', size = SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize = SMALL_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize = MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize = SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize = BIGGER_SIZE)  # fontsize of the figure title

""" transpiration plot """
dt_ = 60 / (24 * 3600)
potential_trans = lambda t, dt: trans * sinusoidal2(t, dt)

# load data
n = len(fnames)
data = [np.load(path + n_ + ".npy")  for n_ in fnames]

fig, ax = plt.subplots(n, 1, figsize = (18, 8))
if n == 1:
    ax = [ax]

for i in range(0, n):
    t = data[i][0]
    y = data[i][1]
    if trans > 0:
        ax[i].plot(t, potential_trans(t, dt_ * np.ones(t.shape)), 'k', label = "potential transpiration")  # potential transpiration
    ax[i].plot(t, y, 'g', label = "actual transpiration ({:s})".format(str_[i]))  # actual transpiration  according to soil model
    ax[i].set_ylabel("transpiration [cm$^3$ day$^{-1}$]")
    ax[i].legend(loc = 'upper left')
    ax2 = ax[i].twinx()
    dt = np.diff(t)
    so = np.array(y)
    cup = np.cumsum(np.multiply(so[:-1], dt))
    ax2.plot(t[1:], cup, 'c--', label = "cumulative transpiration")  # cumulative transpiration (neumann)
    ax2.set_ylabel("cumulative [cm$^3$]")
    ax2.legend(loc = 'center right')
    print(str_[i], "cumulative uptake", cup[-1], "cm3")

ax[i].set_xlabel("Time [d]")

plt.tight_layout()
plt.show()
