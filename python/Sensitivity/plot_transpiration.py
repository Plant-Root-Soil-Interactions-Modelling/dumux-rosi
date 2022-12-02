"""
transpiration plot (one column, number of rows as number of filenames)
"""
import sys; sys.path.append("../../../CPlantBox/");  sys.path.append("../../../CPlantBox/src/python_modules")

import numpy as np
import matplotlib.pyplot as plt

from xylem_flux import sinusoidal2
import evapotranspiration as evap

Kc_maize = 1.2
Kc_soybean = 1.15

# name = "soybean"
# str_ = ["_sra0"]
# sim_time = 87.5
# area = 38 * 5
# range_ = ['1995-03-15 00:00:00', '1995-06-12 11:00:00']
# potential_trans = evap.get_transpiration_beers('data/95.pkl', range_, area, 87.5, evap.lai_soybean, Kc_soybean)
# trans = 1

# name = "test1.0"
# str_ = [""]
# sim_time = 87.5
# potential_trans = lambda t, dt: trans * sinusoidal2(t, dt) * t / sim_time  # soybean
# area = 38 * 5
# range_ = ['1995-03-15 00:00:00', '1995-06-12 11:00:00']
# potential_trans = evap.get_transpiration_beers('data/95.pkl', range_, area, 87.5, evap.lai_soybean, Kc_soybean)
# trans = 1

# name = "maize"
# str_ = ["_cyl0"]
# sim_time = 0.5 * 95
# area = 75 * 15  # cm2
# range_ = ['1995-03-15 00:00:00', '1995-06-20 23:00:00']
# potential_trans = evap.get_transpiration_beers('data/95.pkl', range_, area, 95, evap.lai_maize, Kc_maize)
# trans = 1

name = "maize"
str_ = ["_sra0"]
sim_time = 0.5 * 95
area = 75 * 15  # cm2
range_ = ['1995-03-15 00:00:00', '1995-06-20 23:00:00']
potential_trans = evap.get_transpiration_beers('data/95.pkl', range_, area, 95, evap.lai_maize, Kc_maize)
trans = 1

rs_age = 1

fnames = np.array(["transpiration_" + name + s for s in str_ ])
fnames2 = np.array(["nitrate_" + name + s for s in str_ ])
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
try:
    data2 = [np.load(path + n_ + ".npy")  for n_ in fnames2]
except:
    data2 = None

fig, ax = plt.subplots(n, 1, figsize = (18, 8))
if n == 1:
    ax = [ax]

for i in range(0, n):
    t = data[i][0]
    y = data[i][1]
    if trans > 0:
        ax[i].plot(t, [ -potential_trans(t[i], dt_) / area for i in range(0, t.shape[0]) ], 'k', label = "potential transpiration")  # potential transpiration
    y = np.maximum(y, 0)
    ax[i].plot(t, y / area, 'g', label = "actual transpiration")  # actual transpiration  according to soil model
    if data2 is not None:
        t = data[i][0]
        c = data2[i]
        ax[i].plot(t[::10], 10*c, 'r:', label = "nitrate uptake 0.1*[g/day]")  # actual transpiration  according to soil model
        dt = np.diff(t[::10])
        cuc = np.cumsum(np.multiply(c[:-1], dt))
        print(str_[i], "cumulative nitrate uptake", cuc[-1] , "g")

    ax[i].set_ylabel("transpiration [cm/day]")
    ax[i].legend(loc = 'upper left')
    ax2 = ax[i].twinx()
    dt = np.diff(t)
    so = np.array(y) / area
    cup = np.cumsum(np.multiply(so[:-1], dt))
    ax2.plot(t[1:], cup, 'c--', label = "cumulative transpiration")  # cumulative transpiration (neumann)
    ax2.set_ylabel("cumulative [cm]")
    if data2 is not None:
        ax2.plot(t[10::10], 10*cuc, 'r--', label = "cumulative nitrate uptake 0.1*[g/day]")  # cumulative transpiration (neumann)

    ax2.legend(loc = 'center right')
    print(str_[i], "cumulative water uptake", cup[-1], "cm3")

ax[i].set_xlabel("Time [d]")

plt.tight_layout()
plt.show()
