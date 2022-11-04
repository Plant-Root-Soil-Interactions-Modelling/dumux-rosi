"""
transpiration plot (one column, number of rows as number of filenames)
"""
import sys; sys.path.append("../../../CPlantBox/");  sys.path.append("../../../CPlantBox/src/python_modules")

import numpy as np
import matplotlib.pyplot as plt

from xylem_flux import sinusoidal2
import evapotranspiration as evap

# name = "soybean"
# str_ = ["sra0"]
# sim_time = 0.25 * 87.5
# potential_trans = lambda t, dt: trans * sinusoidal2(t, dt) * t / sim_time  # soybean
# area = 38 * 5
# range_ = ['1995-03-15 00:00:00', '1995-06-10 11:00:00']  # TODO fix range for maize
# potential_trans = evap.get_transpiration('data/95.pkl', range_, area)
# trans = 1

# name = "maize"
# str_ = ["cyl0"]
# sim_time = 0.25 * 95
# area = 75 * 15  # cm2
# range_ = ['1995-03-15 00:00:00', '1995-06-10 11:00:00']  # TODO fix range for maize
# potential_trans = evap.get_transpiration('data/95.pkl', range_, area)
# trans = 1

name = "maize"
str_ = ["sra0"]
sim_time = 0.25 * 95
area = 75 * 15  # cm2
range_ = ['1995-03-15 00:00:00', '1995-06-10 11:00:00']  # TODO fix range for maize
potential_trans = evap.get_transpiration('data/95.pkl', range_, area)
trans = 1

fnames = np.array(["transpiration_" + name + "_" + s for s in str_ ])
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
dt_ = 360 / (24 * 3600)

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
        ax[i].plot(t, [ -potential_trans(t[i], dt_) / area for i in range(0, t.shape[0]) ], 'k', label = "potential transpiration")  # potential transpiration
    ax[i].plot(t, y / area, 'g', label = "actual transpiration ({:s})".format(str_[i]))  # actual transpiration  according to soil model
    ax[i].set_ylabel("transpiration [cm day$^{-1}$]")
    ax[i].legend(loc = 'upper left')
    ax2 = ax[i].twinx()
    dt = np.diff(t)
    so = np.array(y) / area
    cup = np.cumsum(np.multiply(so[:-1], dt))
    ax2.plot(t[1:], cup, 'c--', label = "cumulative transpiration")  # cumulative transpiration (neumann)
    ax2.set_ylabel("cumulative [cm]")
    ax2.legend(loc = 'center right')
    print(str_[i], "cumulative uptake", cup[-1], "cm3")

ax[i].set_xlabel("Time [d]")

plt.tight_layout()
plt.show()
