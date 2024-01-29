"""
     Histograms of Hsr over some points in time 
"""
import matplotlib.pyplot as plt
import numpy as np

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


def plot_hsr(ax, method, dim, plant, soil, outer_method, label_ = "3D", plot_times = [0., 2., 4., 6], ls = [ "-", ":", "-.", "-", ":", "-."]):

    # check with scenario_setup
    l = 150  # cm soil depth
    dx = 1  # cm resolution
    ylim = 120
    days = 14.5

    fnames_hsr = []
    # fnames_hx = []
    for i in range(0, len(outer_method)):
        name = method[i] + "_" + plant + "_" + dim + "_" + soil + "_" + outer_method[i]
        fnames_hsr.append("results/hsr_" + method + "_" + plant + "_" + dim + "_" + soil + "_" + outer_method[i])
        # fnames_hx.append("results/hx_" + method + "_" + plant + "_" + dim + "_" + soil + "_" + outer_method[i])

    cmap = plt.get_cmap('Set1')
    col = cmap([1, 0, 4, 3, 2, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15])  # adjust colors to jans plot

    """ load data """
    n = len(fnames_hsr)
    data = [np.load(n_ + ".npy")  for n_ in fnames_hsr]

    """ sink plot """
    ax[0].set_title("Noon")
    ax[1].set_title("Night")

    """ noon """
    all_min = 1.e6
    all_max = -1.e6
    for i in range(0, n):  # number of files
        hsr_ = data[i]
        peak_id = np.round(hsr_.shape[0] / days * np.array([0.5 + j for j in plot_times]))  # hsr_.shape[0] is the number of saved time steps
        peak_id = peak_id.astype(int)
        for ind, j in enumerate(peak_id):
            all_min = min(all_min, np.min(hsr_[j,:]))
            all_min = max(all_min, -15000)
            all_max = max(all_max, np.max(hsr_[j,:]))
    print("all_min", all_min, "all_max", all_max)

    for i in range(0, n):  # number of files
        hsr_ = data[i]
        peak_id = np.round(hsr_.shape[0] / days * np.array([0.5 + j for j in plot_times]))  # hsr_.shape[0] is the number of saved time steps
        peak_id = peak_id.astype(int)
        for ind, j in enumerate(peak_id):
            lstr = "day {:g} ({:s})".format(plot_times[ind], label_[i])  # method[i]
            # lstr = "{:g}d ({:s})".format(plot_times[ind], outer_method[i])
            # ax[0].plot(sink_[j,:] / cell_volume, soil_z_, label = lstr, color = col[ind], linestyle = ls[i])
            print("Hsr at noon are ranging from ", np.min(hsr_[j,:]), "to ", np.max(hsr_[j,:]), "[cm]")
            ax[0].hist(hsr_[j,:], range = (all_min, all_max), bins = 40, rwidth = 0.9,
                       label = outer_method[i] + ", " + dim + ", day {:g}".format(plot_times[ind]), alpha = 0.5)

        ax[0].legend()

    """ midnight """
    all_min = 1.e6
    all_max = -1.e6
    for i in range(0, n):  # number of files
        hsr_ = data[i]
        redistribution_id = np.round(hsr_.shape[0] / days * np.array([j for j in plot_times]))
        redistribution_id = redistribution_id.astype(int)
        for ind, j in enumerate(redistribution_id):
            all_min = min(all_min, np.min(hsr_[j,:]))
            all_max = max(all_max, np.max(hsr_[j,:]))
    print("all_min", all_min, "all_max", all_max)

    for i in range(0, n):
        hsr_ = data[i]
        redistribution_id = np.round(hsr_.shape[0] / days * np.array([j for j in plot_times]))
        redistribution_id = redistribution_id.astype(int)
        for ind, j in enumerate(redistribution_id):
            # lstr = "{:g}d ({:s})".format(plot_times[ind], outer_method[i])
            lstr = "day {:g} ({:s})".format(plot_times[ind], label_[i])  # method[i], int(ind == 0) +
            # ax[1].plot(sink_[j,:] / cell_volume, soil_z_, label = lstr, color = col[ind], linestyle = ls[i])
            print("Hsr at night are ranging from ", np.min(hsr_[j,:]), "to ", np.max(hsr_[j,:]), "[cm]")
            ax[1].hist(hsr_[j,:], range = (all_min, all_max), bins = 40, rwidth = 0.9,
                       label = outer_method[i] + ", " + dim + ", day {:g}".format(plot_times[ind]), alpha = 0.5)
            # ax[1].hist(hx_[j,:], range = (-1000, 0), bins = 40, rwidth = 0.9, label = ["hx 1D"], alpha = 0.5)

        ax[1].legend()


if __name__ == "__main__":

    # """ maize h_sr over time """
    # fig, ax = plt.subplots(1, 2, figsize = (18, 10))
    # method = "sra"
    # plant = "springbarley"
    # dim = "1D"
    # soil = "hydrus_loam"
    # outer_method = ["voronoi"]
    # plot_hsr(ax, method, dim, plant, soil, outer_method, label_ = dim, plot_times = [1., 4., 6.])  # , 4., 6.
    # plt.savefig('hsr_' + plant + "_" + soil + dim + "_temporal.png")
    # plt.show()

    """ maize h_sr over time """
    fig, ax = plt.subplots(1, 2, figsize = (18, 10))
    method = "sra"
    plant = "maize"
    dim = "1D"
    soil = "hydrus_loam"
    outer_method = ["voronoi", "surface"]
    plot_hsr(ax, method, dim, plant, soil, outer_method, label_ = dim, plot_times = [13.])  # , 4., 6.
    plt.savefig('hsr_' + plant + "_" + soil + dim + ".png")
    plt.show()

    plt.show()

