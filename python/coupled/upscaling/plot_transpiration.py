"""
transpiration plot (one column, number of rows as number of filenames)
"""
import sys; sys.path.append("../../../../CPlantBox");  sys.path.append("../../../../CPlantBox/src")

import numpy as np
import matplotlib.pyplot as plt
from functional.xylem_flux import sinusoidal2

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


def plot_potential(ax, method, plant, dim, soil, outer_method):

    fnames = []
    for i in range(0, len(plant)):
        name = method[i] + "_" + plant[i] + "_" + dim[i] + "_" + soil[i] + "_" + outer_method[i]
        fnames.append("results/transpiration_" + method[i] + "_" + plant[i] + "_" + dim[i] + "_" + soil[i] + "_" + outer_method[i])
    n = len(fnames)
    data = [np.load(n_ + ".npy")  for n_ in fnames]

    potential_trans = lambda t: sinusoidal2(t, dt_)

    # check with scenario_setup
    trans = []
    for p in plant:
        if p == "soybean":
                trans.append(0.5 * 3 * 76)  # potential transpiration cm3/day
        elif p == "maize":
                trans.append(0.5 * 16 * 76)  # potential transpiration cm3/day
        else:
            raise
    dt_ = 360. / (24 * 3600)  # days

    for i in range(0, n):
        t = data[i][0]
        if trans[i] > 0:
            trans_ = trans[i] * potential_trans(t)
            ax[i].plot(t, trans_, 'k', label = "potential transpiration")  # potential transpiration


def plot_transpiration_rows(ax, method, plant, dim, soil, outer_method, ls, label):

    fnames = []
    for i in range(0, len(plant)):
        name = method[i] + "_" + plant[i] + "_" + dim[i] + "_" + soil[i] + "_" + outer_method[i]
        fnames.append("results/transpiration_" + method[i] + "_" + plant[i] + "_" + dim[i] + "_" + soil[i] + "_" + outer_method[i])

    titles = [n[22:] for n in fnames]  # skip 'resuts/transpiration_'

    n = len(fnames)
    data = [np.load(n_ + ".npy")  for n_ in fnames]

    for i in range(0, n):
        t = data[i][0]
        y = data[i][1]
        ax[i].plot(t, y, 'g' + ls, label = " actual transpiration " + label)  # actual transpiration  according to soil model
        ax[i].set_title(titles[i])
        ax[i].set_ylabel("[cm$^3$ day$^{-1}$]")
        ax[i].legend(loc = 'upper left')
        ax2 = ax[i].twinx()
        dt = np.diff(t)
        so = np.array(y)
        cup = np.cumsum(np.multiply(so[:-1], dt))
        ax2.plot(t[1:], cup, 'c--', label = "cumulative transpiration")  # cumulative transpiration (neumann)
        ax2.set_ylabel("cumulative [cm$^3$]")
        ax2.legend(loc = 'center right')
        if plant[i] == "soybean":
            yrange = [0, 1000]  # e.g. for soybean
        else:
            yrange = [0, 4500]  # e.g. for maize
        ax2.set_ylim(yrange)
        print(i, "cumulative uptake", cup[-1], "cm3")

    ax[n - 1].set_xlabel("Time [d]")


if __name__ == "__main__":

    # """ Maize sra vs sraOld """
    # fig, ax = plt.subplots(3, 1, figsize = (12, 14))
    # method = ["sra"] * 3
    # plant = ["maize"] * 3
    # dim = ["1D"] * 3
    # soil = ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]
    # outer_method = ["voronoi"] * 3
    # plot_potential(ax, method, plant, dim, soil, outer_method)
    # plot_transpiration_rows(ax, method, plant, dim, soil, outer_method, "", "(sra)")
    # method = ["sraOld"] * 3
    # plot_transpiration_rows(ax, method, plant, dim, soil, outer_method, ":", "(sraOld)")
    # plt.tight_layout()
    # plt.show()
    # print()
    #
    # """ Soybean sra vs sraOld """
    # fig, ax = plt.subplots(3, 1, figsize = (12, 14))
    # method = ["sra"] * 3
    # plant = ["maize"] * 3
    # dim = ["1D"] * 3
    # soil = ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]
    # outer_method = ["voronoi"] * 3
    # plot_potential(ax, method, plant, dim, soil, outer_method)
    # plot_transpiration_rows(ax, method, plant, dim, soil, outer_method, "", "(sra)")
    # method = ["sraOld"] * 3
    # plot_transpiration_rows(ax, method, plant, dim, soil, outer_method, ":", "(sraOld)")
    # plt.tight_layout()
    # plt.show()
    # print()

    """ Soybean 3D vs 1D """
    fig, ax = plt.subplots(3, 1, figsize = (12, 14))
    method = ["sra"] * 3
    plant = ["maize"] * 3
    dim = ["3D"] * 3
    soil = ["hydrus_loam"] * 3  # , "hydrus_clay", "hydrus_sandyloam"]
    outer_method = ["surface"] * 3
    plot_potential(ax, method, plant, dim, soil, outer_method)
    plot_transpiration_rows(ax, method, plant, dim, soil, outer_method, "", "(3D)")
    dim = ["1D"] * 3
    plot_transpiration_rows(ax, method, plant, dim, soil, outer_method, ":", "(1D)")
    plt.tight_layout()
    plt.show()
    print()

