"""
    Figure 11 and 12         References AAA vs AAB


    Sink plot (noon and midnight), of a 3d soil
    
    modify __main__ to select simualtion result 
"""
import sys; sys.path.append("../../../../CPlantBox");  sys.path.append("../../../../CPlantBox/src")
sys.path.append("../../modules"); sys.path.append("../../../build-cmake/cpp/python_binding/");

import matplotlib.pyplot as plt
import numpy as np

import scenario_setup as setup

from plot_sink1d import *

SMALL_SIZE = 16 + 4
MEDIUM_SIZE = 16 + 4
BIGGER_SIZE = 16 + 8
plt.rc('font', size = SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize = BIGGER_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize = SMALL_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize = BIGGER_SIZE)  # legend fontsize
plt.rc('figure', titlesize = BIGGER_SIZE)  # fontsize of the figure title


def plot_sink3d(ax, method, plant, soil, outer_method, label_ = "3D", plot_times = [0., 2., 4., 6., 13.], ls = ["-", "--", "-.", ":"]):

    soil_names = { "hydrus_loam":"Loam", "hydrus_clay":"Clay", "hydrus_sandyloam": "Sandy loam"}
    dim = ["3D"] * 3

    # if label_ == "BBB":
    #     print("\n 2D \n")
    #     dim = ["2D"] * 3

    # check with scenario_setup
    l = 150  # cm soil depth
    dx = 1  # cm resolution
    ylim = 110
    days = 14.5

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
            min_b, max_b, cell_number = setup.maize_("3D")  ################################## TODO
            layer_volume = 16 * 76 * dx  # cm3
        elif p == "springbarley":
            min_b, max_b, cell_number = setup.springbarley_("3D")
            layer_volume = 3 * 13 * dx  # cm3
        else:
            raise

    # print(cell_number)

    """ load data """
    n = len(fnames)
    data = [np.load(n_ + ".npy")  for n_ in fnames]

    """ sink plot """
    ax[0].set_title(soil_names[soil[0]])
    ax[0].set_ylabel("depth [cm]")
    ax[0].set_xlabel("sink term at noon [1/day]")
    ax[1].set_xlabel("sink term at night [1/day]")
    ax[0].plot([0, 0], [-l, 0.], "k:")
    ax[1].plot([0, 0], [-l, 0.], "k:")

    """ noon """
    for i in range(0, n):

        sink_ = data[i]
        # print("sink_", sink_.shape)
        sink_ = sink_.reshape((sink_.shape[0], cell_number[-1], cell_number[-2], cell_number[-3]))
        # print("sink_", sink_.shape)
        sink_ = np.sum(sink_, 2)
        # print("sink_", sink_.shape)
        sink_ = np.sum(sink_, 2)
        # print("sink_", sink_.shape)

        soil_z_ = np.linspace(-l + dx / 2., -dx / 2., sink_.shape[1])  # segment mids

        peak_id = np.round(sink_.shape[0] / days * np.array([0.5 + i for i in plot_times]))
        peak_id = peak_id.astype(int)

        for ind, j in enumerate(peak_id):
            lstr = "day {:g} ({:s})".format(plot_times[ind], label_[i])  # method[i]
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

        redistribution_id = np.round(sink_.shape[0] / days * np.array([i for i in plot_times]))
        redistribution_id = redistribution_id.astype(int)
        # print("redistribution_id", redistribution_id)

        for ind, j in enumerate(redistribution_id):
            lstr = "day {:g} ({:s})".format(plot_times[ind], label_[i])  # method[i],int(ind == 0) +
            ax[1].plot(sink_[j,:] / layer_volume, soil_z_, label = lstr, color = col[ind], linestyle = ls[i])

        ax[1].set_ylim([-ylim, 0.])
        ax[1].legend()


if __name__ == "__main__":

    """
    References AAA vs AAB
    """
    for s in ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]:  # Figure 11
        fig, ax = plt.subplots(2, 1, figsize = (10, 18))
        method = ["sra"]  # Axx
        plant = ["springbarley"]
        soil = [s]
        outer_method = ["voronoi"]  # xAx
        plot_sink3d(ax, method, plant, soil, outer_method, ["AAA"], [0, 2, 6, 13])
        plot_sink1d(ax, method, plant, soil, outer_method, ["AAB"], [0, 2, 6, 13], ls = ["-."])
        ax[0].set_xlim([0., 0.1])
        ax[1].set_xlim([-0.012, 0.0075])
        ax[1].set_xticks([-0.012, -0.008, -0.004, 0., 0.003, 0.007])
        plt.tight_layout()
        plt.savefig('sink_AAx_springbarley_' + s[7:] + '.png')
    for s in ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]:  # Figure 12
        fig, ax = plt.subplots(2, 1, figsize = (10, 18))
        method = ["sra"]  # Axx
        plant = ["maize"]
        soil = [s]
        outer_method = ["voronoi"]  # xAx
        plot_sink3d(ax, method, plant, soil, outer_method, ["AAA"], [0, 2, 6, 13])
        plot_sink1d(ax, method, plant, soil, outer_method, ["AAB"], [0, 2, 6, 13], ls = ["--"])
        ax[0].set_xlim([0., 0.058])
        ax[1].set_xlim([-0.006, 0.006])
        plt.tight_layout()
        plt.savefig('sink_AAx_maize_' + s[7:] + '.png')

    plt.show()
