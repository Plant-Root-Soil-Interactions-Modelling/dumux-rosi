"""
Sink plot (noon and midnight), of a 1d soil 
"""
import matplotlib.pyplot as plt
import numpy as np

SMALL_SIZE = 16 + 4
MEDIUM_SIZE = 16 + 4
BIGGER_SIZE = 16 + 4
plt.rc('font', size = SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize = SMALL_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize = MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize = SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize = BIGGER_SIZE)  # fontsize of the figure title


def plot_sink_diff(ax, method0, method1, plant, soil, outer_method, plot_times = [0., 2., 4., 6, 13], ls = [ "-", ":", "-.", "-", ":", "-."]):

    soil_names = { "hydrus_loam":"Loam", "hydrus_clay":"Clay", "hydrus_sandyloam": "Sandy loam"}
    dim = "1D"  # DO NOT CHANGE TO 3D, use script plot_sink3d

    # check with scenario_setup
    l = 150  # cm soil depth
    dx = 1  # cm resolution
    days = 14.5

    fnames = []
    name = method0 + "_" + plant + "_" + dim + "_" + soil + "_" + outer_method
    print("Difference plot", name)
    fnames.append("results/sink_" + method0 + "_" + plant + "_" + dim + "_" + soil + "_" + outer_method)
    fnames.append("results/sink_" + method1 + "_" + plant + "_" + dim + "_" + soil + "_" + outer_method)

    cmap = plt.get_cmap('Set1')
    col = cmap([1, 0, 4, 3, 2, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15])  # adjust colors to jans plot

    trans = []
    if plant == "soybean":
            cell_volume = 3 * 76 * dx  # cm3
    elif plant == "maize":
            cell_volume = 16 * 76 * dx  # cm3
    elif plant == "springbarley":
            cell_volume = 3 * 13 * dx  # cm3
    else:
        raise

    """ load data """
    n = len(fnames)
    data = [np.load(n_ + ".npy")  for n_ in fnames]
    sink_diff_ = data[0] - data[1]
    print(sink_diff_.shape[1])
    soil_z_ = np.linspace(-l + dx / 2., -dx / 2., int(sink_diff_.shape[1] / 5))  # segment mids

    """ sink plot """
    ax[0].set_title(soil_names[soil])
    ax[0].set_ylabel("Depth [cm]")
    ax[0].set_xlabel("Difference at noon [1/day]")
    ax[1].set_xlabel("Difference at night [1/day]")
    # ax[0].plot([0, 0], [-l, 0.], "k:")
    # ax[1].plot([0, 0], [-l, 0.], "k:")

    """ noon """
    peak_id = np.round(sink_diff_.shape[0] / days * np.array([0.5 + j for j in plot_times]))
    peak_id = peak_id.astype(int)
    for ind, j in enumerate(peak_id):
        lstr = "day {:g}".format(plot_times[ind])
        width = 5. / len(plot_times)
        diff = np.array([np.mean(sink_diff_[j, k:k + 5]) for k in range(0, sink_diff_.shape[1], 5)])
        ax[0].barh(soil_z_ + np.ones(soil_z_.shape) * ind * width, diff, width, align = 'center', label = lstr)
    ax[0].legend()

    """ midnight """
    redistribution_id = np.round(sink_diff_.shape[0] / days * np.array([j for j in plot_times]))
    redistribution_id = redistribution_id.astype(int)
    for ind, j in enumerate(redistribution_id):
        lstr = "day {:g}".format(plot_times[ind])
        width = 5. / len(plot_times)
        diff = np.array([np.mean(sink_diff_[j, k:k + 5]) for k in range(0, sink_diff_.shape[1], 5)])
        ax[1].barh(soil_z_ + np.ones(soil_z_.shape) * ind * width, diff, width, align = 'center', label = lstr)
    ax[1].legend()


if __name__ == "__main__":

    """ ABB vs. BBB (1d soil, aggregated)"""
    for s in ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]:
        fig, ax = plt.subplots(2, 1, figsize = (10, 18))
        method0 = "sra"
        method1 = "agg"
        plant = "springbarley"
        soil = s
        outer_method = "length"
        plot_sink_diff(ax, method0, method1, plant, soil, outer_method, [0, 2, 6, 13])  #
        ax[0].set_ylim([-110, 0.])
        ax[1].set_ylim([-110, 0.])
        plt.tight_layout()
        plt.savefig('sink_CBB_springbarley_' + s[7:] + '.png')
    for s in ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]:
        fig, ax = plt.subplots(2, 1, figsize = (10, 18))
        method0 = "sra"
        method1 = "agg"
        plant = "maize"
        soil = s
        outer_method = "length"
        plot_sink_diff(ax, method0, method1, plant, soil, outer_method, [0, 2, 6, 13])  #
        ax[0].set_ylim([-95, 0.])
        ax[1].set_ylim([-95, 0.])
        plt.tight_layout()
        plt.savefig('sink_BBB_maize_' + s[7:] + '_diff.png')

    # """ ABB vs. CBB; Parallel 1d vs 1d (plot) """
    # for s in ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]:
    #     fig, ax = plt.subplots(2, 1, figsize = (10, 18))
    #     method0 = "sra"
    #     method1 = "par"
    #     plant = "springbarley"
    #     soil = s
    #     outer_method = "length"
    #     plot_sink_diff(ax, method0, method1, plant, soil, outer_method, [0, 2, 6, 13])
    #     ax[0].set_ylim([-110, 0.])
    #     ax[1].set_ylim([-110, 0.])
    #     plt.tight_layout()
    #     plt.savefig('sink_CBB_springbarley_' + s[7:] + '.png')
    # for s in ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]:
    #     fig, ax = plt.subplots(2, 1, figsize = (10, 18))
    #     method0 = "sra"
    #     method1 = "par"
    #     plant = "maize"
    #     soil = s
    #     outer_method = "length"
    #     plot_sink_diff(ax, method0, method1, plant, soil, outer_method, [0, 2, 6, 13])
    #     ax[0].set_ylim([-95, 0.])
    #     ax[1].set_ylim([-95, 0.])
    #     plt.tight_layout()
    #     plt.savefig('sink_CBB_maize_' + s[7:] + '_diff.png')

    plt.show()

