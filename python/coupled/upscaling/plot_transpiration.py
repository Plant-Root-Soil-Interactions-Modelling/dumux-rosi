"""
transpiration plot (one column, number of rows as number of filenames)
"""
import sys; sys.path.append("../../../../CPlantBox");  sys.path.append("../../../../CPlantBox/src")

import numpy as np
import matplotlib.pyplot as plt
from functional.xylem_flux import sinusoidal2

method = ["sra", "sra", "sra"]
plant = ["maize"] * 3
dim = ["1D"] * 3
soil = ["hydrus_loam"] * 3
outer_method = ["surface"] * 3

fnames = []
for i in range(0, len(plant)):
    name = method[i] + "_" + plant[i] + "_" + dim[i] + "_" + soil[i] + "_" + outer_method[i]
    fnames.append("results/transpiration_" + method[i] + "_" + plant[i] + "_" + dim[i] + "_" + soil[i] + "_" + outer_method[i])

titles = [n[22:] for n in fnames]  # skip 'resuts/transpiration_'

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

# check with scenario_setup
trans = []
for p in plant:
    if p == "soybean":
            trans.append(0.5 * 3 * 75)  # potential transpiration cm3/day
    elif p == "maize":
            trans.append(0.5 * 16 * 75)  # potential transpiration cm3/day
    else:
        raise

dt_ = 360. / (24 * 3600)  # days

""" transpiration plot """
potential_trans = lambda t: sinusoidal2(t, dt_)

# load data
n = len(fnames)
data = [np.load(n_ + ".npy")  for n_ in fnames]

fig, ax = plt.subplots(n, 1, figsize = (14, 12))

for i in range(0, n):
    t = data[i][0]
    y = data[i][1]
    if trans[i] > 0:
        trans_ = trans[i] * potential_trans(t)
        ax[i].plot(t, trans_, 'k', label = "potential transpiration")  # potential transpiration
    ax[i].plot(t, y, 'g', label = " actual transpiration")  # actual transpiration  according to soil model
    ax[i].set_title(titles[i])
    ax[i].set_ylabel("[cm$^3$ day$^{-1}$]")
    # ax[i].set_ylim([0., 1150.])
    ax[i].legend(loc = 'upper left')
    ax2 = ax[i].twinx()
    dt = np.diff(t)
    so = np.array(y)
    cup = np.cumsum(np.multiply(so[:-1], dt))
    ax2.plot(t[1:], cup, 'c--', label = "cumulative transpiration")  # cumulative transpiration (neumann)
    ax2.set_ylabel("cumulative [cm$^3$]")
    # ax2.set_ylim(yrange)
    ax2.legend(loc = 'center right')
    print(i, "cumulative uptake", cup[-1], "cm3")

ax[n - 1].set_xlabel("Time [d]")
plt.tight_layout()
plt.show()
