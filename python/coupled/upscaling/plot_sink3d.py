"""
Sink plot (noon and midnight), of a 3d soil 
"""
import sys; sys.path.append("../../../../CPlantBox");  sys.path.append("../../../../CPlantBox/src")
sys.path.append("../../modules"); sys.path.append("../../../build-cmake/cpp/python_binding/");

import matplotlib.pyplot as plt
import numpy as np

import scenario_setup as setup

from plot_sink1d import *

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


def plot_sink3d(ax, method, plant, soil, outer_method):

    dim = ["3D"] * 3

    # check with scenario_setup
    l = 150  # cm soil depth
    dx = 1  # cm resolution
    ylim = 110
    days = 7.5

    fnames = []
    for i in range(0, len(plant)):
        name = method[i] + "_" + plant[i] + "_" + dim[i] + "_" + soil[i] + "_" + outer_method[i]
        fnames.append("results/sink_" + method[i] + "_" + plant[i] + "_" + dim[i] + "_" + soil[i] + "_" + outer_method[i])

    cmap = plt.get_cmap('Set1')
    col = cmap([1, 0, 4, 3, 2, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15])  # adjust colors to jans plot

    trans = []
    for p in plant:
        if p == "soybean":
            min_b, max_b, cell_number = setup.soybean_("3D")
            layer_volume = 3 * 76 * dx  # cm3
        elif p == "maize":
            min_b, max_b, cell_number = setup.maize_("3D")
            layer_volume = 16 * 76 * dx  # cm3
        elif p == "springbarley":
            min_b, max_b, cell_number = setup.springbarley_("3D")
            layer_volume = 3 * 13 * dx  # cm3
        else:
            raise
    plot_times = [0., 2., 4., 6]

    print(cell_number)

    """ load data """
    n = len(fnames)
    data = [np.load(n_ + ".npy")  for n_ in fnames]

    """ sink plot """
    ax[0].set_title(soil[0][7:])
    ax[0].set_ylabel("depth [cm]")
    ax[0].set_xlabel("sink term at noon [1/day]")
    ax[1].set_xlabel("sink term at night [1/day]")
    ax[0].plot([0, 0], [-l, 0.], "k:")
    ax[1].plot([0, 0], [-l, 0.], "k:")
    ls = ["-", "--", "-.", ":"]

    """ noon """
    for i in range(0, n):

        sink_ = data[i]
        print("sink_", sink_.shape)
        sink_ = sink_.reshape((sink_.shape[0], cell_number[-1], cell_number[-2], cell_number[-3]))
        print("sink_", sink_.shape)
        sink_ = np.sum(sink_, 2)
        print("sink_", sink_.shape)
        sink_ = np.sum(sink_, 2)
        print("sink_", sink_.shape)

        soil_z_ = np.linspace(-l + dx / 2., -dx / 2., sink_.shape[1])  # segment mids

        peak_id = np.round(sink_.shape[0] / days * np.array([0.5 + i for i in plot_times]))
        peak_id = peak_id.astype(int)

        for ind, j in enumerate(peak_id):
            if True:  # ind == 0 or ind == 3:
                lstr = "{:g}d ({:s})".format(plot_times[ind], "3D")  # method[i]
                ax[0].plot(sink_[j,:] / layer_volume, soil_z_, label = lstr, color = col[ind], linestyle = ls[i])

        ax[0].set_ylim([-ylim, 0.])
        ax[0].legend()

    """ midnight """
    for i in range(0, n):

        sink_ = data[i]
        sink_ = sink_.reshape((sink_.shape[0], cell_number[-1], cell_number[-2], cell_number[-3]))
        sink_ = np.sum(sink_, 2)
        sink_ = np.sum(sink_, 2)

        soil_z_ = np.linspace(-l + dx / 2., -dx / 2., sink_.shape[1])  # segment mids

        redistribution_id = np.round(sink_.shape[0] / days * np.array([1 + i for i in plot_times]))
        redistribution_id = redistribution_id.astype(int)
        print("redistribution_id", redistribution_id)

        for ind, j in enumerate(redistribution_id):
            if True:  # ind == 0 or ind == 3:
                lstr = "{:g}d ({:s})".format(plot_times[ind], "3D")  # method[i],int(ind == 0) +
                ax[1].plot(sink_[j,:] / layer_volume, soil_z_, label = lstr, color = col[ind], linestyle = ls[i])

        ax[1].set_ylim([-ylim, 0.])
        ax[1].legend()


if __name__ == "__main__":

    """ springbarley """  # check L92, 114 (should be True:)
    for s in ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]:
        fig, ax = plt.subplots(2, 1, figsize = (10, 18))
        method = ["agg"]
        plant = ["springbarley"]
        soil = [s]
        outer_method = ["voronoi"]
        plot_sink3d(ax, method, plant, soil, outer_method)
        ax[0].set_xlim([0., 0.115])
        ax[1].set_xlim([-0.012, 0.005])
        # plot_sink1d(ax, method, plant, soil, outer_method)
        plt.tight_layout()
        plt.savefig('sink3d_springbarley_' + s[7:] + '.png')

    # """ springbarley 1d vs 3d"""  # check L92, 114 (should be only 2 time points)
    # for s in ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]:
    #     fig, ax = plt.subplots(2, 1, figsize = (10, 18))
    #     method = ["sra"]
    #     plant = ["springbarley"]
    #     soil = [s]
    #     outer_method = ["surface"]
    #     plot_sink3d(ax, method, plant, soil, outer_method)
    #     ax[0].set_xlim([0., 0.115])
    #     ax[1].set_xlim([-0.012, 0.005])
    #     plot_sink1d(ax, method, plant, soil, outer_method)
    #     plt.tight_layout()
    #     plt.savefig('sink1d_springbarley_' + s[7:] + '.png')

    # """ maize 3d """
    # for s in ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]:
    #     fig, ax = plt.subplots(2, 1, figsize = (10, 18))
    #     method = ["sra"]
    #     plant = ["maize"]
    #     soil = [s]
    #     outer_method = ["surface"]
    #     plot_sink3d(ax, method, plant, soil, outer_method)
    #     ax[0].set_xlim([0., 0.038])
    #     ax[1].set_xlim([-0.0025, 0.0023])
    #     plt.tight_layout()
    #     plt.savefig('sink3d_maize_' + s[7:] + '.png')

    # """ maize 1d vs 3d"""
    # for s in ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]:
    #     fig, ax = plt.subplots(2, 1, figsize = (10, 18))
    #     method = ["sra"]
    #     plant = ["maize"]
    #     soil = [s]
    #     outer_method = ["surface"]
    #     plot_sink3d(ax, method, plant, soil, outer_method)
    #     plot_sink1d(ax, method, plant, soil, outer_method)
    #     ax[0].set_xlim([0., 0.038])
    #     ax[1].set_xlim([-0.0025, 0.0023])
    #     plt.tight_layout()
    #     plt.savefig('sink1d_maize_' + s[7:] + '.png')

    plt.show()

    # method = ["sra", "sraOld", "sra"]
    # plant = ["maize"] * 3
    # dim = ["1D"] * 3
    # soil = ["hydrus_clay", "hydrus_loam", "hydrus_sandyloam"]
    # outer_method = ["surface"] * 3
