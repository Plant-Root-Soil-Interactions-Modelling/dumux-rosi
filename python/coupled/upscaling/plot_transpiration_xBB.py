"""
    transpiration plot (one column, number of rows as number of filenames)
    
    modify __main__ to select simualtion result    
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

    cup_ = []
    cup2_ = []
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
        ax2[i].plot(t[1:], cup, 'c' + ls, label = "cumulative transpiration " + label)  # cumulative transpiration (neumann)
        ax2[i].set_ylabel("cumulative [cm$^3$]")
        ax2[i].legend(loc = 'lower right')
        print(i, "cumulative uptake " + label, cup[-1], "cm3")
        cup_.append(cup[-1])
        cup2_.append(cup[cup.shape[0] // 2])
        # ax[i].set_xlim([0, 7])  # data ranges until 7.5
    ax[n - 1].set_xlabel("Time [d]")

    return np.array(cup_), np.array(cup2_)


if __name__ == "__main__":

    # """
    # Aggregated model
    # """
    # print("\nAggregated 3d vs. Reference 3d (plot)")
    # fig, ax = plt.subplots(3, 1, figsize = (12, 14))
    # ax2 = [ ax[i].twinx() for i in range(0, len(ax)) ]
    # method = ["sra"] * 3
    # plant = ["springbarley"] * 3
    # dim = ["3D"] * 3
    # soil = ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]
    # outer_method = ["voronoi"] * 3
    # plot_potential(ax, method, plant, dim, soil, outer_method)
    # cup_ref, cup2_ref = plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, ":", "(ref)")
    # method = ["agg"] * 3
    # cup_, cup2_ = plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, "", "(agg)")
    # plt.savefig('transpiration_agg3d_springbarley.png')
    # print("\nSpringbarley: percental error in cumulative uptake comapared to reference solution")
    # print("1 week ", 100.*(np.ones(np.shape(cup2_)) - np.divide(cup2_, cup2_ref)), "%")
    # print("2 weeks", 100.*(np.ones(np.shape(cup_)) - np.divide(cup_, cup_ref)), "%")
    # ##
    # print()
    # fig, ax = plt.subplots(3, 1, figsize = (12, 14))
    # ax2 = [ ax[i].twinx() for i in range(0, len(ax)) ]
    # method = ["sra"] * 3
    # plant = ["maize"] * 3
    # dim = ["3D"] * 3
    # soil = ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]
    # outer_method = ["voronoi"] * 3
    # plot_potential(ax, method, plant, dim, soil, outer_method)
    # cup_ref, cup2_ref = plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, ":", "(ref)")
    # method = ["agg"] * 3
    # cup_, cup2_ = plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, "", "(agg)")
    # plt.savefig('transpiration_agg3d_maize.png')
    # print("\nMaize: percental error in cumulative uptake comapared to reference solution")
    # print("1 week ", 100.*(np.ones(np.shape(cup2_)) - np.divide(cup2_, cup2_ref)), "%")
    # print("2 weeks", 100.*(np.ones(np.shape(cup_)) - np.divide(cup_, cup_ref)), "%")
    # #
    # print("\nAggregated 1d vs. Reference 3d ")
    # fig, ax = plt.subplots(3, 1, figsize = (12, 14))
    # ax2 = [ ax[i].twinx() for i in range(0, len(ax)) ]
    # method = ["sra"] * 3
    # plant = ["springbarley"] * 3
    # dim = ["3D"] * 3
    # soil = ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]
    # outer_method = ["voronoi"] * 3
    # plot_potential(ax, method, plant, dim, soil, outer_method)
    # cup_ref, cup2_ref = plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, ":", "(ref)")
    # dim = ["1D"] * 3
    # method = ["agg"] * 3
    # cup_, cup2_ = plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, "", "(agg 1D)")
    # plt.savefig('transpiration_agg1d_springbarley.png')
    # print("\nSpringbarley: percental error in cumulative uptake comapared to reference soluWer kümmert sich bei einem neuen Scrum-Team um die Ressourcen?tion")
    # print("1 week ", 100.*(np.ones(np.shape(cup2_)) - np.divide(cup2_, cup2_ref)), "%")
    # print("2 weeks", 100.*(np.ones(np.shape(cup_)) - np.divide(cup_, cup_ref)), "%")
    # # ##
    # print()
    # fig, ax = plt.subplots(3, 1, figsize = (12, 14))
    # ax2 = [ ax[i].twinx() for i in range(0, len(ax)) ]
    # method = ["sra"] * 3
    # plant = ["maize"] * 3
    # dim = ["3D"] * 3
    # soil = ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]
    # outer_method = ["voronoi"] * 3
    # plot_potential(ax, method, plant, dim, soil, outer_method)
    # cup_ref, cup2_ref = plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, ":", "(ref)")
    # dim = ["1D"] * 3
    # method = ["agg"] * 3
    # cup_ = plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, "", "(agg 1D)")
    # plt.savefig('transpiration_agg1d_maize.png')
    # print("Maize: percental error in cumulative uptake comapared to reference solution")
    # print("1 week ", 100.*(np.ones(np.shape(cup2_)) - np.divide(cup2_, cup2_ref)), "%")
    # print("2 weeks", 100.*(np.ones(np.shape(cup_)) - np.divide(cup_, cup_ref)), "%")

    # print("\nABB vs. BBB")
    # fig, ax = plt.subplots(3, 1, figsize = (12, 14))
    # ax2 = [ ax[i].twinx() for i in range(0, len(ax)) ]
    # method = ["sra"] * 3
    # plant = ["springbarley"] * 3
    # dim = ["1D"] * 3
    # soil = ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]
    # outer_method = ["length"] * 3
    # plot_potential(ax, method, plant, dim, soil, outer_method)
    # cup_ref, cup2_ref = plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, "", "(ABB)")
    # method = ["agg"] * 3
    # cup_, cup2_ = plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, ":", "(BBB)")
    # plt.savefig('transpiration_BBB_springbarley.png')
    # print("\nSpringbarley: percental error in cumulative uptake comapared to reference solution")
    # print("1 week ", 100.*(np.ones(np.shape(cup2_)) - np.divide(cup2_, cup2_ref)), "%")
    # print("2 weeks", 100.*(np.ones(np.shape(cup_)) - np.divide(cup_, cup_ref)), "%")
    # # ##
    # print()
    # fig, ax = plt.subplots(3, 1, figsize = (12, 14))
    # ax2 = [ ax[i].twinx() for i in range(0, len(ax)) ]
    # method = ["sra"] * 3
    # plant = ["maize"] * 3
    # dim = ["1D"] * 3
    # soil = ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]
    # outer_method = ["length"] * 3
    # plot_potential(ax, method, plant, dim, soil, outer_method)
    # cup_ref, cup2_ref = plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, "", "(ABB)")
    # method = ["agg"] * 3
    # cup_, cup2_ = plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, ":", "(BBB)")
    # plt.savefig('transpiration_BBB_maize.png')
    # print("\nMaize: percental error in cumulative uptake comapared to reference solution")
    # print("1 week ", 100.*(np.ones(np.shape(cup2_)) - np.divide(cup2_, cup2_ref)), "%")
    # print("2 weeks", 100.*(np.ones(np.shape(cup_)) - np.divide(cup_, cup_ref)), "%")

    """
    Parallel root system
    """
    # print("\nParallel 3d vs. Reference 3d ")
    # fig, ax = plt.subplots(3, 1, figsize = (12, 14))
    # ax2 = [ ax[i].twinx() for i in range(0, len(ax)) ]
    # method = ["sra"] * 3
    # plant = ["springbarley"] * 3
    # dim = ["3D"] * 3
    # soil = ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]
    # outer_method = ["surface"] * 3  ##########################################################################################
    # plot_potential(ax, method, plant, dim, soil, outer_method)
    # cup_ref, cup2_ref = plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, ":", "(ref)")
    # method = ["par"] * 3
    # cup_, cup2_ = plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, "", "(par)")
    # plt.savefig('transpiration_par3d_springbarley.png')
    # print("\nSpringbarley: percental error in cumulative uptake comapared to reference solution")
    # print("1 week ", 100.*(np.ones(np.shape(cup2_)) - np.divide(cup2_, cup2_ref)), "%")
    # print("2 weeks", 100.*(np.ones(np.shape(cup_)) - np.divide(cup_, cup_ref)), "%")
    # # ##
    # print()
    # fig, ax = plt.subplots(3, 1, figsize = (12, 14))
    # ax2 = [ ax[i].twinx() for i in range(0, len(ax)) ]
    # method = ["sra"] * 3
    # plant = ["maize"] * 3
    # dim = ["3D"] * 3
    # soil = ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]
    # outer_method = ["surface"] * 3  ##########################################################################################
    # plot_potential(ax, method, plant, dim, soil, outer_method)
    # cup_ref, cup2_ref = plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, ":", "(ref)")
    # method = ["par"] * 3
    # cup_, cup2_ = plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, "", "(par)")
    # plt.savefig('transpiration_par3d_maize.png')
    # print("\nMaize: percental error in cumulative uptake comapared to reference solution")
    # print("1 week ", 100.*(np.ones(np.shape(cup2_)) - np.divide(cup2_, cup2_ref)), "%")
    # print("2 weeks", 100.*(np.ones(np.shape(cup_)) - np.divide(cup_, cup_ref)), "%")
    # plt.tight_layout()
    # plt.show()
    # print()

    # print("\nParallel 1d vs. Reference 3d ")
    # fig, ax = plt.subplots(3, 1, figsize = (12, 14))
    # ax2 = [ ax[i].twinx() for i in range(0, len(ax)) ]
    # method = ["sra"] * 3
    # plant = ["springbarley"] * 3
    # dim = ["3D"] * 3
    # soil = ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]
    # outer_method = ["voronoi"] * 3
    # plot_potential(ax, method, plant, dim, soil, outer_method)
    # cup_ref, cup2_ref = plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, ":", "(ref)")
    # dim = ["1D"] * 3
    # method = ["par"] * 3
    # cup_, cup2_ = plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, "", "(par 1D)")
    # plt.savefig('transpiration_par1d_springbarley.png')
    # print("\nSpringbarley: percental error in cumulative uptake comapared to reference soluWer kümmert sich bei einem neuen Scrum-Team um die Ressourcen?tion")
    # print("1 week ", 100.*(np.ones(np.shape(cup2_)) - np.divide(cup2_, cup2_ref)), "%")
    # print("2 weeks", 100.*(np.ones(np.shape(cup_)) - np.divide(cup_, cup_ref)), "%")
    # # ##
    # print()
    # fig, ax = plt.subplots(3, 1, figsize = (12, 14))
    # ax2 = [ ax[i].twinx() for i in range(0, len(ax)) ]
    # method = ["sra"] * 3
    # plant = ["maize"] * 3
    # dim = ["3D"] * 3
    # soil = ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]
    # outer_method = ["length"] * 3
    # plot_potential(ax, method, plant, dim, soil, outer_method)
    # cup_ref, cup2_ref = plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, ":", "(ref)")
    # dim = ["1D"] * 3
    # method = ["par"] * 3
    # cup_, cup2_ = plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, "", "(par 1D)")
    # plt.savefig('transpiration_par1d_maize.png')
    # print("Maize: percental error in cumulative uptake comapared to reference solution")
    # print("1 week ", 100.*(np.ones(np.shape(cup2_)) - np.divide(cup2_, cup2_ref)), "%")
    # print("2 weeks", 100.*(np.ones(np.shape(cup_)) - np.divide(cup_, cup_ref)), "%")
    # plt.tight_layout()
    # plt.show()
    # print()

    print("\nABB vs. CBB")
    fig, ax = plt.subplots(3, 1, figsize = (12, 14))
    ax2 = [ ax[i].twinx() for i in range(0, len(ax)) ]
    method = ["sra"] * 3
    plant = ["springbarley"] * 3
    dim = ["1D"] * 3
    soil = ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]
    outer_method = ["length"] * 3
    plot_potential(ax, method, plant, dim, soil, outer_method)
    cup_ref, cup2_ref = plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, "", "(ABB)")
    method = ["par"] * 3
    cup_, cup2_ = plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, ":", "(CBB)")
    plt.savefig('transpiration_CBB_springbarley.png')
    print("\nSpringbarley: percental error in cumulative uptake comapared to reference solution")
    print("1 week ", 100.*(np.ones(np.shape(cup2_)) - np.divide(cup2_, cup2_ref)), "%")
    print("2 weeks", 100.*(np.ones(np.shape(cup_)) - np.divide(cup_, cup_ref)), "%")
    # ##
    print()
    fig, ax = plt.subplots(3, 1, figsize = (12, 14))
    ax2 = [ ax[i].twinx() for i in range(0, len(ax)) ]
    method = ["sra"] * 3
    plant = ["maize"] * 3
    dim = ["1D"] * 3
    soil = ["hydrus_loam", "hydrus_clay", "hydrus_sandyloam"]
    outer_method = ["length"] * 3
    plot_potential(ax, method, plant, dim, soil, outer_method)
    cup_ref, cup2_ref = plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, "", "(ABB)")
    method = ["par"] * 3
    cup_, cup2_ = plot_transpiration_rows(ax, ax2, method, plant, dim, soil, outer_method, ":", "(CBB)")
    plt.savefig('transpiration_CBB_maize.png')
    print("\nMaize: percental error in cumulative uptake comapared to reference solution")
    print("1 week ", 100.*(np.ones(np.shape(cup2_)) - np.divide(cup2_, cup2_ref)), "%")
    print("2 weeks", 100.*(np.ones(np.shape(cup_)) - np.divide(cup_, cup_ref)), "%")

    plt.tight_layout()
    plt.show()