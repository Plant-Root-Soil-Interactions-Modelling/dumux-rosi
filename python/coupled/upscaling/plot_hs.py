"""
    Sink plot (noon and midnight), of a 3d soil
    
    modify __main__ to select simualtion result 
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


def plot_hs(ax, method, plant, soil, outer_method, label_ = "3D", plot_times = [0., 2., 6., 13.], ls = ["-", "--", "-.", ":"]):  #  4., 6., 13.

    dim = ["3D"] * 3

    # check with scenario_setup
    l = 150  # cm soil depth
    dx = 1  # cm resolution
    ylim = 110
    days = 14.5

    fnames, fnames1D = [], []
    for i in range(0, len(plant)):
        fnames.append("results/hs_" + method[i] + "_" + plant[i] + "_3D_" + soil[i] + "_" + outer_method[i])
        fnames1D.append("results/hs_" + method[i] + "_" + plant[i] + "_1D_" + soil[i] + "_" + outer_method[i])

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
            print("plot_hs() unknonw plant name")
            raise

    # print(cell_number)

    """ load data """
    n = len(fnames)
    data = [np.load(n_ + ".npy")  for n_ in fnames]
    data1D = [np.load(n_ + ".npy")  for n_ in fnames1D]

    """ sink plot """
    ax[0].set_title(soil[0][7:])
    ax[0].set_ylabel("depth [cm]")
    ax[0].set_xlabel("soil matric potential [cm]")

    """ noon """
    for i in range(0, n):

        hs_ = data[i]
        hs_ = hs_.reshape((hs_.shape[0], cell_number[-1], cell_number[-2], cell_number[-3]))
        hs2_ = np.mean(hs_, 2)
        hs2_ = np.mean(hs2_, 2)
        hs1D = data1D[i]

        soil_z_ = np.linspace(-l + dx / 2., -dx / 2., hs_.shape[1])  # segment mids

        peak_id = np.round(hs_.shape[0] / days * np.array([i for i in plot_times]))  # 0.5 +
        peak_id = peak_id.astype(int)

        for ind, j in enumerate(peak_id):
            lstr = "day {:g} ({:s})".format(plot_times[ind], label_[i])  # method[i]
            ax[0].plot(hs2_[j,:], soil_z_, label = "3D, " + lstr, color = col[ind])
            # ax[0].plot(hs1D[j,:], soil_z_, label = "1D, " + lstr, color = col[ind], linestyle = "--")

            for k in range(0, len(soil_z_)):
                d_ = hs_[j, k,:].flat
                ax[0].plot(d_, [soil_z_[k]] * len(d_), '*', color = col[ind], alpha = 0.1 * 0.1)
        # ax[0].set_ylim([-ylim, 0.])
        ax[0].legend()


if __name__ == "__main__":

    for s in ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]:
        fig, ax = plt.subplots(1, 1, figsize = (10, 18))
        ax = [ax]
        method = ["sra"]
        plant = ["maize"]
        soil = [s]
        outer_method = ["voronoi"]
        plot_hs(ax, method, plant, soil, outer_method, ["ref"])
        # ax[0].set_xlim([0., 0.038])
        # ax[1].set_xlim([-0.0025, 0.0025])
        plt.tight_layout()
        plt.savefig('hs_3d_' + plant[0] + "_" + outer_method[0] + '_' + s[7:] + '.png')

    """ 
    References 3d, 1d 
    """
    # """ Reference 3d (plot) """
    # for s in ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]:
    #     fig, ax = plt.subplots(2, 1, figsize = (10, 18))
    #     method = ["sra"]
    #     plant = ["springbarley"]
    #     soil = [s]
    #     outer_method = ["voronoi"]
    #     plot_sink3d(ax, method, plant, soil, outer_method, ["ref"])
    #     ax[0].set_xlim([0., 0.12])
    #     ax[1].set_xlim([-0.012, 0.0075])
    #     plt.tight_layout()
    #     plt.savefig('sink_3d_springbarley_' + s[7:] + '.png')
    # for s in ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]:
    #     fig, ax = plt.subplots(2, 1, figsize = (10, 18))
    #     method = ["sra"]
    #     plant = ["maize"]
    #     soil = [s]
    #     outer_method = ["voronoi"]
    #     plot_sink3d(ax, method, plant, soil, outer_method, ["ref"])
    #     ax[0].set_xlim([0., 0.038])
    #     ax[1].set_xlim([-0.0025, 0.0025])
    #     plt.tight_layout()
    #     plt.savefig('sink_3d_maize_' + s[7:] + '.png')
    #
    # """ Reference 1d vs 3d (plot) """
    # for s in ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]:
    #     fig, ax = plt.subplots(2, 1, figsize = (10, 18))
    #     method = ["sra"]
    #     plant = ["springbarley"]
    #     soil = [s]
    #     outer_method = ["voronoi"]
    #     plot_sink1d(ax, method, plant, soil, outer_method, ["1D"], [0, 2, 6, 13])
    #     plot_sink3d(ax, method, plant, soil, outer_method, ["ref"], [0, 2, 6, 13], ls = ["--"])
    #     ax[0].set_xlim([0., 0.12])
    #     ax[1].set_xlim([-0.012, 0.0075])
    #     plt.tight_layout()
    #     plt.savefig('sink_1d_springbarley_' + s[7:] + '.png')
    # for s in ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]:
    #     fig, ax = plt.subplots(2, 1, figsize = (10, 18))
    #     method = ["sra"]
    #     plant = ["maize"]
    #     soil = [s]
    #     outer_method = ["voronoi"]
    #     plot_sink1d(ax, method, plant, soil, outer_method, ["1D"], [0, 2, 6, 13])
    #     plot_sink3d(ax, method, plant, soil, outer_method, ["ref"], [0, 2, 6, 13], ls = ["--"])
    #     ax[0].set_xlim([0., 0.058])
    #     ax[1].set_xlim([-0.006, 0.006])
    #     plt.tight_layout()
    #     plt.savefig('sink_1d_maize_' + s[7:] + '.png')
    #
    # """
    # Aggregated model
    # """
    # """  Aggregaed 3d vs Reference 3d (plot) """
    # for s in ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]:
    #     fig, ax = plt.subplots(2, 1, figsize = (10, 18))
    #     method = ["agg", "sra"]
    #     plant = ["springbarley"] * 2
    #     soil = [s] * 2
    #     outer_method = ["voronoi"] * 2
    #     plot_sink3d(ax, method, plant, soil, outer_method, ["agg", "ref"], [0, 2, 6, 13])
    #     ax[0].set_xlim([0., 0.12])
    #     ax[1].set_xlim([-0.012, 0.0075])
    #     plt.tight_layout()
    #     plt.savefig('sink_agg3d_springbarley_' + s[7:] + '.png')
    # for s in ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]:
    #     fig, ax = plt.subplots(2, 1, figsize = (10, 18))
    #     method = ["agg", "sra"]
    #     plant = ["maize"] * 2
    #     soil = [s] * 2
    #     outer_method = ["voronoi"] * 2
    #     plot_sink3d(ax, method, plant, soil, outer_method, ["agg", "ref"], [0, 2, 6, 13])
    #     ax[0].set_xlim([0., 0.038])
    #     ax[1].set_xlim([-0.0025, 0.0025])
    #     plt.tight_layout()
    #     plt.savefig('sink_agg3d_maize_' + s[7:] + '.png')
    #
    # """ Aggregatd 1d vs 3d (plot) """
    # for s in ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]:
    #     fig, ax = plt.subplots(2, 1, figsize = (10, 18))
    #     method = ["agg"]
    #     plant = ["springbarley"]
    #     soil = [s]
    #     outer_method = ["voronoi"]
    #     plot_sink1d(ax, method, plant, soil, outer_method, ["agg 1D"], [0, 2, 6, 13])
    #     method = ["sra"]
    #     plot_sink3d(ax, method, plant, soil, outer_method, ["ref"], [0, 2, 6, 13], ls = ["--"])
    #     ax[0].set_xlim([0., 0.12])
    #     ax[1].set_xlim([-0.012, 0.0075])
    #     plt.tight_layout()
    #     plt.savefig('sink_agg1d_springbarley_' + s[7:] + '.png')
    # for s in ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]:
    #     fig, ax = plt.subplots(2, 1, figsize = (10, 18))
    #     method = ["agg"]
    #     plant = ["maize"]
    #     soil = [s]
    #     outer_method = ["voronoi"]
    #     plot_sink1d(ax, method, plant, soil, outer_method, ["agg 1D"], [0, 2, 6, 13])
    #     method = ["sra"]
    #     plot_sink3d(ax, method, plant, soil, outer_method, ["ref"], [0, 2, 6, 13], ls = ["--"])
    #     ax[0].set_xlim([0., 0.058])
    #     ax[1].set_xlim([-0.006, 0.006])
    #     plt.tight_layout()
    #     plt.savefig('sink_agg1d_maize_' + s[7:] + '.png')
    #
    # """ Aggregatd 1d vs 1d (plot) """
    # for s in ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]:
    #     fig, ax = plt.subplots(2, 1, figsize = (10, 18))
    #     method = ["agg"]
    #     plant = ["springbarley"]
    #     soil = [s]
    #     outer_method = ["voronoi"]
    #     plot_sink1d(ax, method, plant, soil, outer_method, ["agg 1D"], [0, 2, 6, 13])
    #     method = ["sra"]
    #     plot_sink1d(ax, method, plant, soil, outer_method, ["ref 1D"], [0, 2, 6, 13], ls = ["--"])
    #     ax[0].set_xlim([0., 0.12])
    #     ax[1].set_xlim([-0.012, 0.0075])
    #     plt.tight_layout()
    #     plt.savefig('sink_agg1d1d_springbarley_' + s[7:] + '.png')
    # for s in ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]:
    #     fig, ax = plt.subplots(2, 1, figsize = (10, 18))
    #     method = ["agg"]
    #     plant = ["maize"]
    #     soil = [s]
    #     outer_method = ["voronoi"]
    #     plot_sink1d(ax, method, plant, soil, outer_method, ["agg 1D"], [0, 2, 6, 13])
    #     method = ["sra"]
    #     plot_sink1d(ax, method, plant, soil, outer_method, ["ref 1D"], [0, 2, 6, 13], ls = ["--"])
    #     ax[0].set_xlim([0., 0.058])
    #     ax[1].set_xlim([-0.006, 0.006])
    #     plt.tight_layout()
    #     plt.savefig('sink_agg1d1d_maize_' + s[7:] + '.png')

    """
    Parallel root system
    """
    """  Parallel 3d vs Reference 3d (plot) """
    # for s in ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]:
    #     fig, ax = plt.subplots(2, 1, figsize = (10, 18))
    #     method = ["par", "sra"]
    #     plant = ["springbarley"] * 2
    #     soil = [s] * 2
    #     outer_method = ["length"] * 2
    #     plot_sink3d(ax, method, plant, soil, outer_method, ["par", "ref"], [0, 2, 6, 13])
    #     ax[0].set_xlim([0., 0.12])
    #     ax[1].set_xlim([-0.012, 0.0075])
    #     plt.tight_layout()
    #     plt.savefig('sink_par3d_springbarley_' + s[7:] + '.png')
    # for s in ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]:
    #     fig, ax = plt.subplots(2, 1, figsize = (10, 18))
    #     method = ["par", "sra"]
    #     plant = ["maize"] * 2
    #     soil = [s] * 2
    #     outer_method = ["surface"] * 2
    #     plot_sink3d(ax, method, plant, soil, outer_method, ["par", "ref"], [0, 2, 6, 13])
    #     ax[0].set_xlim([0., 0.038])
    #     ax[1].set_xlim([-0.0025, 0.0025])
    #     plt.tight_layout()
    #     plt.savefig('sink_par3d_maize_' + s[7:] + '.png')
    # """ Parallel 1d vs 3d (plot) """
    # for s in ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]:
    #     fig, ax = plt.subplots(2, 1, figsize = (10, 18))
    #     method = ["par"]
    #     plant = ["springbarley"]
    #     soil = [s]
    #     outer_method = ["voronoi"]
    #     plot_sink1d(ax, method, plant, soil, outer_method, ["par 1D"], [0, 2, 6, 13])
    #     method = ["sra"]
    #     plot_sink3d(ax, method, plant, soil, outer_method, ["ref"], [0, 2, 6, 13], ls = ["--"])
    #     ax[0].set_xlim([0., 0.12])
    #     ax[1].set_xlim([-0.012, 0.0075])
    #     plt.tight_layout()
    #     plt.savefig('sink_par1d_springbarley_' + s[7:] + '.png')
    # for s in ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]:
    #     fig, ax = plt.subplots(2, 1, figsize = (10, 18))
    #     method = ["par"]
    #     plant = ["maize"]
    #     soil = [s]
    #     outer_method = ["voronoi"]
    #     plot_sink1d(ax, method, plant, soil, outer_method, ["par 1D"], [0, 2, 6, 13])
    #     method = ["sra"]
    #     plot_sink3d(ax, method, plant, soil, outer_method, ["ref"], [0, 2, 6, 13], ls = ["--"])
    #     ax[0].set_xlim([0., 0.058])
    #     ax[1].set_xlim([-0.006, 0.006])
    #     plt.tight_layout()
    #     plt.savefig('sink_par1d_maize_' + s[7:] + '.png')
    """ Parallel 1d vs 1d (plot) """
    # for s in ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]:
    #     fig, ax = plt.subplots(2, 1, figsize = (10, 18))
    #     method = ["par"]
    #     plant = ["springbarley"]
    #     soil = [s]
    #     outer_method = ["voronoi"]
    #     plot_sink1d(ax, method, plant, soil, outer_method, ["par 1D"], [0, 2, 6, 13])
    #     method = ["sra"]
    #     plot_sink1d(ax, method, plant, soil, outer_method, ["ref 1D"], [0, 2, 6, 13], ls = ["--"])
    #     ax[0].set_xlim([0., 0.12])
    #     ax[1].set_xlim([-0.012, 0.0075])
    #     plt.tight_layout()
    #     plt.savefig('sink_par1d1d_springbarley_' + s[7:] + '.png')
    # for s in ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]:
    #     fig, ax = plt.subplots(2, 1, figsize = (10, 18))
    #     method = ["par"]
    #     plant = ["maize"]
    #     soil = [s]
    #     outer_method = ["voronoi"]
    #     plot_sink1d(ax, method, plant, soil, outer_method, ["par 1D"], [0, 2, 6, 13])
    #     method = ["sra"]
    #     plot_sink1d(ax, method, plant, soil, outer_method, ["ref 1D"], [0, 2, 6, 13], ls = ["--"])
    #     ax[0].set_xlim([0., 0.058])
    #     ax[1].set_xlim([-0.006, 0.006])
    #     plt.tight_layout()
    #     plt.savefig('sink_par1d1d_maize_' + s[7:] + '.png')

    plt.show()
