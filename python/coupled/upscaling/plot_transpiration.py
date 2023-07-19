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
        elif p == "springbarley":
                trans.append(0.5 * 3 * 13)  # potential transpiration cm3/day
        else:
            raise
    dt_ = 360. / (24 * 3600)  # days

    for i in range(0, n):
        t = data[i][0]
        if trans[i] > 0:
            trans_ = trans[i] * potential_trans(t)
            ax[i].plot(t, trans_, 'k', label = "potential transpiration")  # potential transpiration


def plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, ls, label):

    soil_ = ["Loam", "Clay", "Sandy loam"]

    fnames = []
    for i in range(0, len(plant)):
        name = method[i] + "_" + plant[i] + "_" + dim[i] + "_" + soil[i] + "_" + outer_method[i]
        fnames.append("results/transpiration_" + method[i] + "_" + plant[i] + "_" + dim[i] + "_" + soil[i] + "_" + outer_method[i])

    titles = [n[22:] for n in fnames]  # skip 'resuts/transpiration_'
    titles = [n[:-8] for n in titles]

    n = len(fnames)
    data = [np.load(n_ + ".npy")  for n_ in fnames]

    for i in range(0, n):
        t = data[i][0]
        y = data[i][2]
        ax[i].plot(t, y, 'g' + ls, label = " actual transpiration " + label)  # actual transpiration  according to soil model
        ax[i].set_title(soil_[i])
        # ax[i].set_title(titles[i])
        ax[i].set_ylabel("[cm$^3$ day$^{-1}$]")
        ax[i].legend(loc = 'upper left')
        dt = np.diff(t)
        so = np.array(y)
        cup = np.cumsum(np.multiply(so[:-1], dt))
        ax2[i].plot(t[1:], cup, 'c' + ls, label = "cumulative transpiration" + label)  # cumulative transpiration (neumann)
        ax2[i].set_ylabel("cumulative [cm$^3$]")
        ax2[i].legend(loc = 'lower right')
        print(i, "cumulative uptake" + label, cup[-1], "cm3")
        ax[i].set_xlim([0, 7])  # data ranges until 7.5

    ax[n - 1].set_xlabel("Time [d]")


if __name__ == "__main__":

    """ Soil types 3d"""

    # fig, ax = plt.subplots(3, 1, figsize = (12, 14))
    # ax2 = [ ax[i].twinx() for i in range(0, len(ax)) ]
    # method = ["sra"] * 3
    # plant = ["maize"] * 3
    # dim = ["3D"] * 3
    # soil = ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]
    # outer_method = ["surface"] * 3
    # plot_potential(ax, method, plant, dim, soil, outer_method)
    # plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, "", "(surf)")
    # outer_method = ["voronoi"] * 3
    # plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, ":", "(vor)")
    # plt.savefig('transpiration3d_maize.png')

    # fig, ax = plt.subplots(3, 1, figsize = (12, 14))
    # ax2 = [ ax[i].twinx() for i in range(0, len(ax)) ]
    # method = ["sra"] * 3
    # plant = ["soybean"] * 3
    # dim = ["3D"] * 3
    # soil = ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]
    # outer_method = ["surface"] * 3
    # plot_potential(ax, method, plant, dim, soil, outer_method)
    # plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, "", "(surf)")
    # outer_method = ["voronoi"] * 3
    # plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, ":", "(vor)")
    # plt.savefig('transpiration3d_soybean.png')

    # fig, ax = plt.subplots(3, 1, figsize = (12, 14))
    # ax2 = [ ax[i].twinx() for i in range(0, len(ax)) ]
    # method = ["sra"] * 3
    # plant = ["springbarley"] * 3
    # dim = ["3D"] * 3
    # soil = ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]
    # outer_method = ["surface"] * 3
    # plot_potential(ax, method, plant, dim, soil, outer_method)
    # plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, "", "(surf)")
    # outer_method = ["voronoi"] * 3
    # plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, ":", "(vor)")
    # plt.savefig('transpiration3d_springbarley.png')

    """ Soil types 1d"""
    fig, ax = plt.subplots(3, 1, figsize = (12, 14))
    ax2 = [ ax[i].twinx() for i in range(0, len(ax)) ]
    method = ["agg"] * 3
    plant = ["springbarley"] * 3
    dim = ["1D"] * 3
    soil = ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]
    outer_method = ["voronoi"] * 3
    plot_potential(ax, method, plant, dim, soil, outer_method)
    plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, "", "(vor)")
    # outer_method = ["voronoi"] * 3
    # plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, ":", "(vor)")
    plt.savefig('transpiration1d_maize.png')

    # fig, ax = plt.subplots(3, 1, figsize = (12, 14))
    # ax2 = [ ax[i].twinx() for i in range(0, len(ax)) ]
    # method = ["sra"] * 3
    # plant = ["soybean"] * 3
    # dim = ["1D"] * 3
    # soil = ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]
    # outer_method = ["surface"] * 3
    # plot_potential(ax, method, plant, dim, soil, outer_method)
    # plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, "", "(surf)")
    # outer_method = ["voronoi"] * 3
    # plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, ":", "(vor)")
    # plt.savefig('transpiration1d_soybean.png')

    # fig, ax = plt.subplots(3, 1, figsize = (12, 14))
    # ax2 = [ ax[i].twinx() for i in range(0, len(ax)) ]
    # method = ["sra"] * 3
    # plant = ["springbarley"] * 3
    # dim = ["1D"] * 3
    # soil = ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]
    # outer_method = ["surface"] * 3
    # plot_potential(ax, method, plant, dim, soil, outer_method)
    # plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, "", "(surf)")
    # outer_method = ["voronoi"] * 3
    # plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, ":", "(vor)")
    # plt.savefig('transpiration1d_springbarley.png')

    # """ Soybean 3D vs 1D """
    # fig, ax = plt.subplots(3, 1, figsize = (12, 14))
    # ax2 = [ ax[i].twinx() for i in range(0, len(ax)) ]
    # method = ["sra"] * 3
    # plant = ["maize"] * 3
    # dim = ["3D"] * 3
    # soil = ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]
    # outer_method = ["surface"] * 3
    # plot_potential(ax, method, plant, dim, soil, outer_method)
    # plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, "", "(3D, surf)")
    # dim = ["1D"] * 3
    # outer_method = ["voronoi"] * 3
    # plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, ":", "(1D, vor)")
    # outer_method = ["surface"] * 3
    # plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, "--", "(1D, surf)")
    # plt.savefig('transpiration3dvs1d_maize.png')

    # fig, ax = plt.subplots(3, 1, figsize = (12, 14))
    # ax2 = [ ax[i].twinx() for i in range(0, len(ax)) ]
    # method = ["sra"] * 3
    # plant = ["soybean"] * 3
    # dim = ["3D"] * 3
    # soil = ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]
    # outer_method = ["surface"] * 3
    # plot_potential(ax, method, plant, dim, soil, outer_method)
    # plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, "", "(3D, surf)")
    # dim = ["1D"] * 3
    # outer_method = ["voronoi"] * 3
    # plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, ":", "(1D, vor)")
    # outer_method = ["surface"] * 3
    # plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, "--", "(1D, surf)")
    # plt.savefig('transpiration3dvs1d_soybean.png')

    # fig, ax = plt.subplots(3, 1, figsize = (12, 14))
    # ax2 = [ ax[i].twinx() for i in range(0, len(ax)) ]
    # method = ["sra"] * 3
    # plant = ["springbarley"] * 3
    # dim = ["3D"] * 3
    # soil = ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]
    # outer_method = ["surface"] * 3
    # plot_potential(ax, method, plant, dim, soil, outer_method)
    # plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, "", "(3D, surf)")
    # dim = ["1D"] * 3
    # # outer_method = ["voronoi"] * 3
    # # plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, ":", "(1D, vor)")
    # outer_method = ["surface"] * 3
    # plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, "--", "(1D, surf)")
    # plt.savefig('transpiration3dvs1d_springbarley.png')

    """ todo: benchmark with old implementation """
    # """ (benchmark) Maize sra vs sraOld """
    # method = ["sra"] * 3
    # plant = ["maize"] * 3
    # dim = ["3D"] * 3
    # soil = ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]
    # outer_method = ["voronoi"] * 3
    # plot_potential(ax, method, plant, dim, soil, outer_method)
    # plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, "", "(sra)")
    # method = ["sraOld"] * 3
    # plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, ":", "(sraOld)")

    # """ (benchmark) Soybean sra vs sraOld """
    # method = ["sra"] * 3
    # plant = ["maize"] * 3
    # dim = ["1D"] * 3
    # soil = ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]
    # outer_method = ["voronoi"] * 3
    # plot_potential(ax, method, plant, dim, soil, outer_method)
    # plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, "", "(sra)")
    # method = ["sraOld"] * 3
    # plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, ":", "(sraOld)")

    plt.tight_layout()
    plt.show()
    print()

