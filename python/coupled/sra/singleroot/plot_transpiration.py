"""
transpiration plot (one column, number of rows as number of filenames)
"""

import numpy as np
import matplotlib.pyplot as plt

add_str = "_dry"  # _wet, _dry

# # Constant scenarios
# fnames = ["results/transpiration_singleroot_cyl_constkrkx" + add_str,
#           "results/transpiration_singleroot_sra_constkrkx" + add_str,
#           "results/transpiration_singleroot_agg_constkrkx" + add_str]
# trans = 0  # no potential transpiration

# Dynamic scanrios
fnames = ["results/transpiration_singleroot_cyl_dynamic_constkrkx" + add_str,
          "results/transpiration_singleroot_sra_dynamic_constkrkx" + add_str,
          "results/transpiration_singleroot_agg_dynamic_constkrkx" + add_str]
trans = 0.5  # potential transpiration cm3/day

titles = ["Rhizosphere", "Steady rate", "Aggregated steady rate"]

SMALL_SIZE = 16
MEDIUM_SIZE = 16
BIGGER_SIZE = 16
plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


def sinusoidal(t):
    return np.sin(2. * np.pi * np.array(t) - 0.5 * np.pi) + 1.


""" transpiration plot """
potential_trans = lambda t: trans * sinusoidal(t)

# load data
n = len(fnames)
data = [np.loadtxt(n_, delimiter=';')  for n_ in fnames]

fig, ax = plt.subplots(n, 1, figsize=(14, 12))

for i in range(0, n):
    t = data[i][0]
    y = data[i][1]
    if trans > 0:
        ax[i].plot(t, potential_trans(t), 'k', label="potential transpiration")  # potential transpiration
    ax[i].plot(t, y, 'g', label=" actual transpiration")  # actual transpiration  according to soil model
    # ax[i].set_xlabel("Time [d]")
    ax[i].set_title(titles[i] + " (" + add_str[1:] + ")")
    ax[i].set_ylabel("[cm$^3$ day$^{-1}$]")
    # ax[i].set_ylim([0., 1150.])
    ax[i].legend(loc='upper left')
    ax2 = ax[i].twinx()
    dt = np.diff(t)
    so = np.array(y)
    cup = np.cumsum(np.multiply(so[:-1], dt))
    ax2.plot(t[1:], cup, 'c--', label="cumulative transpiration")  # cumulative transpiration (neumann)
    ax2.set_ylabel("cumulative [cm$^3$]")
    # ax2.set_ylim(yrange)
    ax2.legend(loc='center right')
    print(i, "cumulative uptake", cup[-1], "cm3")

plt.tight_layout()
plt.show()