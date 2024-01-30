"""
     Histograms of Hsr over some points in time, sorted into bins 
"""
import sys; sys.path.append("../../../build-cmake/cpp/python_binding/"); sys.path.append("../../../python/modules/");
sys.path.append("../../../../CPlantBox"); sys.path.append("../../../../CPlantBox/src")

import plot_rootsystem as pr

from functional.Perirhizal import *

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


def to_zrange(x:list, other:list = [], min_:float = -np.inf, max_:float = np.inf):
    """ returns ths list with x values within min_ and max_ (and drops nans) 
        drops the same entries from @param other list
    """
    y, z = [], []
    for i, x_ in enumerate(x):
        if (not np.isnan(x_)) and x_ >= min_ and x_ <= max_:
            y.append(x_)
            if other:
                z.append(other[i])
    return y, z


def plot_hsr_bins(ax, method, dim, plant, soil, outer_method, label_ = "3D", plot_times = [0., 2., 4., 6], ls = [ "-", ":", "-.", "-", ":", "-."]):

    outer_r_bins = False  # False: bins per detph
    vis_hsr = False  # False: Hs; True: Hsr

    # check with scenario_setup
    l = 150  # cm soil depth
    dx = 1  # cm resolution
    ylim = 120
    days = 14.5

    # get outer radii & segment z coordinates
    outer_ = []
    z_ = []
    for method_ in outer_method:
        if outer_r_bins == True:
            outer_.append(pr.get_outer_radius(plant, dim, method_))

        r, mapping, length, a, surf, z = pr.get_rootsystem(plant, dim)
        z_.append(z)
        min_z_ = np.min(z_)

    """ load data """
    fnames_hsr = []
    for i in range(0, len(outer_method)):
        name = method[i] + "_" + plant + "_" + dim + "_" + soil + "_" + outer_method[i]
        if vis_hsr:
            fnames_hsr.append("results/hsr_" + method + "_" + plant + "_" + dim + "_" + soil + "_" + outer_method[i])
        else:
            fnames_hsr.append("results/hs_" + method + "_" + plant + "_" + dim + "_" + soil + "_" + outer_method[i])

    n = len(fnames_hsr)
    data = [np.load(n_ + ".npy")  for n_ in fnames_hsr]
    if not vis_hsr:
        data2 = []
        for i, d in enumerate(data):
            d2 = []
            print(d.shape)
            for j in range(0, d.shape[0]):
                d2.append(r.get_hs(list(data[i][j,:])))  # converts from hs per cell to hs per segment
            data2.append(np.array(d2))
    data = data2  # rename

    cmap = plt.get_cmap('Set1')
    col = cmap([1, 0, 4, 3, 2, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15])  # adjust colors to jans plot
    noon = True

    """ noon or night """
    bin_range = 2
    bin_size = bin_range / 10
    bin_size_z = min_z_ / 10
    all_min = [1.e6] * 10
    all_max = [-1.e6] * 10
    if outer_r_bins:
        for i in range(0, n):  # number of files
            hsr_ = data[i]
            if not vis_hsr:
                hsr_ = r.get_hs(list(hsr_))  # converts from hs per cell to hs per segment
            peak_id = np.round(hsr_.shape[0] / days * np.array([0.5 * noon + j for j in plot_times]))  # hsr_.shape[0] is the number of saved time steps
            peak_id = peak_id.astype(int)
            for ind, j in enumerate(peak_id):
                for i_ in range(0, 10):
                    min_or = i * bin_size
                    max_or = (i + 1) * bin_size
                    # print(len(outer_[i]), len(list(hsr_[j,:])))
                    outer_bin, hsr_bin = PerirhizalPython.to_range_(None, outer_[i], list(hsr_[j,:]), min_or, max_or)
                    if hsr_bin:
                        all_min[i_] = min(all_min[i_], np.min(hsr_bin))
                        all_min[i_] = max(all_min[i_], -15000)
                        all_max[i_] = max(all_max[i_], np.max(hsr_bin))
                    else:
                        print("bin", min_or, max_or, "has no entries")
            print("all_min", all_min, "all_max", all_max)

    for i in range(0, n):  # number of files

        hsr_ = data[i]
        peak_id = np.round(hsr_.shape[0] / days * np.array([0.5 * noon + j for j in plot_times]))  # hsr_.shape[0] is the number of saved time steps
        peak_id = peak_id.astype(int)

        for ind, j in enumerate(peak_id):

            for i_ in range(0, 10):  # bins

                if outer_r_bins:
                    min_or = i_ * bin_size
                    max_or = (i_ + 1) * bin_size
                    outer_bin, hsr_bin = PerirhizalPython.to_range_(None, outer_[i], list(hsr_[j,:]), min_or, max_or)
                    print(i, i_, i_ % 5, i_ // 5, "hsr_bin", min_or, max_or, len(hsr_bin), "ranging", all_min[i_], all_max[i_])
                else:
                    max_z = i_ * bin_size_z
                    min_z = (i_ + 1) * bin_size_z
                    # print(len(z_[i]))
                    # print(len(list(hsr_[j,:])))
                    outer_bin, hsr_bin = to_zrange(z_[i], list(hsr_[j,:]), min_z, max_z)
                    print(i, i_, i_ % 5, i_ // 5, "hsr_bin", min_z, max_z, len(hsr_bin), "ranging", min_z_, 0)

                lstr = "{:s}, {:s} , day {:g}".format(outer_method[i], dim, plot_times[ind])
                if outer_r_bins:
                    ax[i_ % 5, i_ // 5].set_title("Bin ({:g}-{:g}) cm outer radius".format(min_or, max_or))
                    if len(hsr_bin) > 0:
                        ax[i_ % 5, i_ // 5].hist(np.array(hsr_bin), range = (all_min[i_], all_max[i_]), bins = 10, rwidth = 0.9, alpha = 0.5, label = lstr)
                else:
                    if len(hsr_bin) > 0:
                        ax[i_ % 5, i_ // 5].hist(np.array(hsr_bin), range = (-3000, 0), bins = 10, rwidth = 0.9, alpha = 0.5, label = lstr)
                    ax[i_ % 5, i_ // 5].set_title("Bin ({:g}-{:g}) cm depth".format(min_z, max_z))

        for i_ in range(0, 10):
            ax[i_ % 5, i_ // 5].legend()


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

    """ maize voronoi vs. uniform """

    fig, ax = plt.subplots(5, 2, figsize = (18, 14))
    method = "sra"
    plant = "maize"
    dim = "1D"
    soil = "hydrus_loam"
    outer_method = ["voronoi", "length"]
    plot_hsr_bins(ax, method, dim, plant, soil, outer_method, label_ = dim, plot_times = [13.])  # , 4., 6.
    plt.tight_layout()
    plt.savefig('hsr_bin_' + plant + "_" + soil + dim + "z.png")
    plt.show()

