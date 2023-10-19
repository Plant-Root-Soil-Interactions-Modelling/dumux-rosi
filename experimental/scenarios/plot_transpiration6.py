"""
transpiration plot, two columns (for dry and wet scenario), number of rows as number of filenames
"""

import numpy as np
import matplotlib.pyplot as plt

names_dry = ["results/transpiration_" + "small_agg" + "_dry",
          "results/transpiration_" + "small_agg" + "_dry"]
names_wet = ["results/transpiration_" + "small_agg" + "_wet",
          "results/transpiration_" + "small_agg" + "_wet"]
titles = ["Aggregated steady rate", "Aggregated steady rate"]  # "Rhizosphere", "Steady rate", "Aggregated steady rate"
trans = 0.5 * 15 * 75  # potential transpiration cm3/day

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


def sinusoidal(t):
    return np.sin(2. * np.pi * np.array(t) - 0.5 * np.pi) + 1.


""" transpiration plot """
potential_trans = lambda t: trans * sinusoidal(t)

# load data
n = len(names_dry)
data_dry = [np.loadtxt(n_, delimiter = ';')  for n_ in names_dry]
data_wet = [np.loadtxt(n_, delimiter = ';')  for n_ in names_wet]

fig, ax = plt.subplots(n, 2, figsize = (28, 12))

for i in range(0, n):
    t = data_wet[i][0]
    y = data_wet[i][1]
    if trans > 0:
        ax[i][0].plot(t, potential_trans(t), 'k', label = "potential transpiration")  # potential transpiration
    ax[i][0].plot(t, y, 'g', label = " actual transpiration")  # actual transpiration  according to soil model
    # ax[i].set_xlabel("Time [d]")
    ax[i][0].set_title(titles[i] + " (wet)")
    ax[i][0].set_ylabel("[cm$^3$ d$^{-1}$]")
    # ax[i].set_ylim([0., 1150.])
    ax[i][0].legend(loc = 'upper left')
    ax2 = ax[i][0].twinx()
    dt = np.diff(t)
    so = np.array(y)
    cup = np.cumsum(np.multiply(so[:-1], dt))
    ax2.plot(t[1:], cup, 'c--', label = "cumulative transpiration")  # cumulative transpiration (neumann)
    ax2.set_ylabel("cumulative [cm$^3$]")
    # ax2.set_ylim(yrange)
    ax2.legend(loc = 'lower right')
    print("wet cumulative uptake", cup[-1], "cm3")

    t = data_dry[i][0]
    y = data_dry[i][1]
    if trans > 0:
        ax[i][1].plot(t, potential_trans(t), 'k', label = "potential transpiration")  # potential transpiration
    ax[i][1].plot(t, y, 'g', label = " actual transpiration")  # actual transpiration  according to soil model
    # ax[i].set_xlabel("Time [d]")
    ax[i][1].set_title(titles[i] + " (dry)")
    ax[i][1].set_ylabel("[cm$^3$ d$^{-1}$]")
    # ax[i].set_ylim([0., 1150.])
    ax[i][1].legend(loc = 'upper left')
    ax2 = ax[i][1].twinx()
    dt = np.diff(t)
    so = np.array(y)
    cup = np.cumsum(np.multiply(so[:-1], dt))
    ax2.plot(t[1:], cup, 'c--', label = "cumulative transpiration")  # cumulative transpiration (neumann)
    ax2.set_ylabel("cumulative [cm$^3$]")
    # ax2.set_ylim(yrange)
    ax2.legend(loc = 'lower right')
    print("dry cumulative uptake", cup[-1], "cm3\n")

plt.tight_layout()
plt.show()
