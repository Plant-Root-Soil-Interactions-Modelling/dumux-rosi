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
    fnames_hx = []
    for i in range(0, len(plant)):
        name = method[i] + "_" + plant[i] + "_" + dim[i] + "_" + soil[i] + "_" + outer_method[i]
        fnames_hsr.append("results/hsr_" + method[i] + "_" + plant[i] + "_" + dim[i] + "_" + soil[i] + "_" + outer_method[i])
        fnames_hx.append("results/hx_" + method[i] + "_" + plant[i] + "_" + dim[i] + "_" + soil[i] + "_" + outer_method[i])

    cmap = plt.get_cmap('Set1')
    col = cmap([1, 0, 4, 3, 2, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15])  # adjust colors to jans plot

    """ load data """
    n = len(fnames_hsr)
    data = [np.load(n_ + ".npy")  for n_ in fnames_hsr]
    data_hx = [np.load(n_ + ".npy")  for n_ in fnames_hx]

    """ sink plot """
    ax[0].set_title("Noon")
    ax[1].set_title("Night")

    """ noon """
    for i in range(0, n):  # number of files

        hsr_ = data[i]
        hx_ = data_hx[i]

        peak_id = np.round(hsr_.shape[0] / days * np.array([0.5 + j for j in plot_times]))  # hsr_.shape[0] is the number of saved time steps
        peak_id = peak_id.astype(int)
        for ind, j in enumerate(peak_id):
            lstr = "day {:g} ({:s})".format(plot_times[ind], label_[i])  # method[i]
            # lstr = "{:g}d ({:s})".format(plot_times[ind], outer_method[i])
            # ax[0].plot(sink_[j,:] / cell_volume, soil_z_, label = lstr, color = col[ind], linestyle = ls[i])
            print("Hsr at noon are ranging from ", np.min(hsr_[j,:]), "to ", np.max(hsr_[j,:]), "[cm]")
            print(len(hsr_[j,:]))
            print(len(hx_[j,:]))
            ax[0].hist(hsr_[j,:], range = (-2000, 0), bins = 40, rwidth = 0.9, label = ["hsr 1D"], alpha = 0.5)
            ax[0].hist(hx_[j,:], range = (-2000, 0), bins = 40, rwidth = 0.9, label = ["hx 1D"], alpha = 0.5)

        ax[0].legend()

    """ midnight """
    for i in range(0, n):

        hsr_ = data[i]

        redistribution_id = np.round(hsr_.shape[0] / days * np.array([j for j in plot_times]))
        redistribution_id = redistribution_id.astype(int)

        for ind, j in enumerate(redistribution_id):
            # lstr = "{:g}d ({:s})".format(plot_times[ind], outer_method[i])
            lstr = "day {:g} ({:s})".format(plot_times[ind], label_[i])  # method[i], int(ind == 0) +
            # ax[1].plot(sink_[j,:] / cell_volume, soil_z_, label = lstr, color = col[ind], linestyle = ls[i])
            print("Hsr at night are ranging from ", np.min(hsr_[j,:]), "to ", np.max(hsr_[j,:]), "[cm]")
            ax[1].hist(hsr_[j,:], range = (-1000, 0), bins = 40, rwidth = 0.9, label = ["hsr 1D"], alpha = 0.5)
            ax[1].hist(hx_[j,:], range = (-1000, 0), bins = 40, rwidth = 0.9, label = ["hx 1D"], alpha = 0.5)

        ax[1].legend()


if __name__ == "__main__":

    """ use plot_sink3d.py, plot_sink1d() is called from there """

    # print("1d vs 3d comparison plots were done with plot_sink3d.py")

    """ springbarley """
    fig, ax = plt.subplots(1, 2, figsize = (18, 10))
    method = ["sra"]
    plant = ["maize"]
    dim = ["1D"]
    soil = ["hydrus_loam"]
    outer_method = ["length"]
    plot_hsr(ax, method, dim, plant, soil, outer_method, label_ = dim, plot_times = [2.])  # , 4., 6.
    plt.savefig('hsr_springbarley_loam.png')
    plt.show()

    # fig, ax = plt.subplots(1, 2, figsize = (18, 10))
    # method = ["agg"]
    # plant = ["springbarley"]
    # soil = ["hydrus_clay"]
    # outer_method = ["voronoi"]  # , "voronoi"]
    # plot_sink1d(ax, method, plant, soil, outer_method)
    # plt.savefig('sink1d_springbarley_clay.png')
    # plt.show()
    #
    # fig, ax = plt.subplots(1, 2, figsize = (18, 10))
    # method = ["agg"]
    # plant = ["springbarley"]
    # soil = ["hydrus_sandyloam"]
    # outer_method = ["voronoi"]  # , "voronoi"]
    # plot_sink1d(ax, method, plant, soil, outer_method)
    # plt.savefig('sink1d_springbarley_sandyloam.png')

    """ maize """
    # fig, ax = plt.subplots(1, 2, figsize = (18, 10))
    # method = ["sra"]
    # plant = ["maize"]
    # soil = ["hydrus_loam"]
    # outer_method = ["surface"]  # , "voronoi"]
    #
    # plot_sink1d(ax, method, plant,  soil,outer_method)
    # plt.savefig('sink1d_maize_loam.png')
    #
    # fig, ax = plt.subplots(1, 2, figsize = (18, 10))
    # method = ["sra"]
    # plant = ["maize"]
    # soil = ["hydrus_clay"]
    # outer_method = ["surface"]  # , "voronoi"]
    # plot_sink1d(ax, method, plant,  soil,outer_method)
    # plt.savefig('sink1d_maize_clay.png')
    #
    # fig, ax = plt.subplots(1, 2, figsize = (18, 10))
    # method = ["sra"]
    # plant = ["maize"]
    # soil = ["hydrus_sandyloam"]
    # outer_method = ["surface"]  # , "voronoi"]
    # plot_sink1d(ax, method, plant,  soil,outer_method)
    # plt.savefig('sink1d_maize_sandyloam.png')
    #

    plt.show()

