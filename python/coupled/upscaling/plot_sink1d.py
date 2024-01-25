"""
Sink plot (noon and midnight), of a 1d soil 
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


def plot_sink1d(ax, method, plant, soil, outer_method, label_ = "3D", plot_times = [0., 2., 4., 6, 13], ls = [ "-", ":", "-.", "-", ":", "-."]):

    dim = ["1D"] * 3  # DO NOT CHANGE TO 3D, use script plot_sink3d

    # check with scenario_setup
    l = 150  # cm soil depth
    dx = 1  # cm resolution
    ylim = 120
    days = 14.5

    fnames = []
    for i in range(0, len(outer_method)):
        name = method[i] + "_" + plant[i] + "_" + dim[i] + "_" + soil[i] + "_" + outer_method[i]
        fnames.append("results/sink_" + method[i] + "_" + plant[i] + "_" + dim[i] + "_" + soil[i] + "_" + outer_method[i])

    cmap = plt.get_cmap('Set1')
    col = cmap([1, 0, 4, 3, 2, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15])  # adjust colors to jans plot

    trans = []
    for p in plant:
        if p == "soybean":
                cell_volume = 3 * 76 * dx  # cm3
        elif p == "maize":
                cell_volume = 16 * 76 * dx  # cm3
        elif p == "springbarley":
                cell_volume = 3 * 13 * dx  # cm3
        else:
            raise

    """ load data """
    n = len(fnames)
    data = [np.load(n_ + ".npy")  for n_ in fnames]

    """ sink plot """
    ax[0].set_title(soil[0])
    ax[0].set_ylabel("depth [cm]")
    ax[0].set_xlabel("sink term at noon [1/day]")
    ax[1].set_xlabel("sink term at night [1/day]")
    ax[0].plot([0, 0], [-l, 0.], "k:")
    ax[1].plot([0, 0], [-l, 0.], "k:")

    """ noon """
    for i in range(0, n):

        sink_ = data[i]
        soil_z_ = np.linspace(-l + dx / 2., -dx / 2., sink_.shape[1])  # segment mids

        peak_id = np.round(sink_.shape[0] / days * np.array([0.5 + j for j in plot_times]))
        peak_id = peak_id.astype(int)

        for ind, j in enumerate(peak_id):
            lstr = "day {:g} ({:s})".format(plot_times[ind], label_[i])  # method[i]
            # lstr = "{:g}d ({:s})".format(plot_times[ind], outer_method[i])
            ax[0].plot(sink_[j,:] / cell_volume, soil_z_, label = lstr, color = col[ind], linestyle = ls[i])

        ax[0].set_ylim([-ylim, 0.])
        ax[0].legend()

    """ midnight """
    for i in range(0, n):

        sink_ = data[i]
        soil_z_ = np.linspace(-l + dx / 2., -dx / 2., sink_.shape[1])  # segment mids

        redistribution_id = np.round(sink_.shape[0] / days * np.array([j for j in plot_times]))
        redistribution_id = redistribution_id.astype(int)

        for ind, j in enumerate(redistribution_id):
            # lstr = "{:g}d ({:s})".format(plot_times[ind], outer_method[i])
            lstr = "day {:g} ({:s})".format(plot_times[ind], label_[i])  # method[i], int(ind == 0) +
            ax[1].plot(sink_[j,:] / cell_volume, soil_z_, label = lstr, color = col[ind], linestyle = ls[i])

        ax[1].set_ylim([-ylim, 0.])
        ax[1].legend()


if __name__ == "__main__":

    """ use plot_sink3d.py, plot_sink1d() is called from there """

    # print("1d vs 3d comparison plots were done with plot_sink3d.py")

    """ springbarley """
    # fig, ax = plt.subplots(1, 2, figsize = (18, 10))
    # method = ["sra", "agg"]
    # plant = ["springbarley"] * 2
    # soil = ["hydrus_loam"] * 2
    # outer_method = ["voronoi"] * 2  # , "voronoi"]
    # plot_sink1d(ax, method, plant, soil, outer_method, label=)
    # plt.savefig('sink1d_springbarley_loam.png')
    # plt.show()

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
    fig, ax = plt.subplots(2, 1, figsize = (10, 18))
    method = ["sra"] * 3
    plant = ["maize"] * 3
    soil = ["hydrus_sandyloam"] * 3
    outer_method = ["voronoi", "length", "surface"]  # , "voronoi"]
    plot_sink1d(ax, method, plant, soil, outer_method, label_ = ["voronoi", "length", "surface"], plot_times = [2., 6., 9., 13])
    plt.savefig('sink1d_' + soil[0] + '.png')
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

