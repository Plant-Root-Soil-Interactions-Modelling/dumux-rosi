""" 
    plots histograms of the outer radius of the perirhizal zone for various discretisations
"""
import sys; sys.path.append("../../build-cmake/cpp/python_binding/"); sys.path.append("../../python/modules/");
sys.path.append("../../../CPlantBox"); sys.path.append("../../../CPlantBox/src")

import plantbox as pb
from functional.Perirhizal import *
import visualisation.vtk_plot as vp
import scenario_setup as scenario

import vtk
import numpy as np
import matplotlib.pyplot as plt

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
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']


def get_outer_radii(rootsystem, type_str, dim):

    r, outer_radii, rs_age, trans, wilting_point, soil, s, sra_table_lookup, mapping = scenario.set_scenario(rootsystem, dim, -200, "hydrus_loam", type_str)

    r = r.rs
    seg2cell = np.array(r.getSegmentMapper())  # r.seg2cell as list

    ana = pb.SegmentAnalyser(r.mappedSegments())
    ana.addData("outer_r", outer_radii)
    ana.addData("length", ana.getParameter("length"))
    # organic_layer = pb.SDF_Cuboid(pb.Vector3d([-1e6, -1e6, -organic]), pb.Vector3d([1e6, 1e6, max_b[2]]))
    topsoil_layer = pb.SDF_Cuboid(pb.Vector3d([-1e6, -1e6, -topsoil]), pb.Vector3d([1e6, 1e6, 0.]))
    subsoil_layer = pb.SDF_Cuboid(pb.Vector3d([-1e6, -1e6, -150.]), pb.Vector3d([1e6, 1e6, -topsoil]))
    # ana0 = pb.SegmentAnalyser(ana)
    # ana0.crop(organic_layer)
    # ana0.pack()
    ana1 = pb.SegmentAnalyser(ana)
    ana1.crop(topsoil_layer)
    ana1.pack()
    ana2 = pb.SegmentAnalyser(ana)
    ana2.crop(subsoil_layer)
    ana2.pack()
    # outer0 = ana0.data["outer_r"]
    outer1 = ana1.data["outer_r"]
    outer2 = ana2.data["outer_r"]
    length1 = ana1.data["length"]
    length2 = ana2.data["length"]
    print("outer1", np.nanmin(outer1), np.nanmax(outer1), "median", np.nanmedian(outer1), "mean_top", np.nanmean(outer1), np.nanstd(outer1))
    print("outer2", np.nanmin(outer2), np.nanmax(outer2), "median", np.nanmedian(outer2), "mean_sub", np.nanmean(outer2), np.nanstd(outer2))  #
    print(len(outer_radii), len(outer1) + len(outer2))  # len(outer0) +

    return outer1, outer2, seg2cell, length1, length2


fig, axes = plt.subplots(2, 3, figsize = (18, 12))

type_str = "length"
topsoil = 30  # 10 * 2.54
subsoil = 150  # 30 * 2.54

""" sring barley """
plant = "springbarley"
titles = []
titles.append("Spring barley 3D grid [1 cm]$^3$")
titles.append("Spring barley 1D grid [1 cm]$^1$")
cell_numbers = []
cell_numbers.append([13, 3, 150])
cell_numbers.append([1, 1, 150])
dim_ = ["3D", "1D"]

stats = []
for i, dim in enumerate(dim_):

    outer1, outer2, seg2cell, length1, length2 = get_outer_radii(plant, type_str, dim)
    stats.append([np.min(outer1), np.max(outer1), np.median(outer1), np.mean(outer1), np.std(outer1),
                  np.min(outer2), np.max(outer2), np.median(outer2), np.mean(outer2), np.std(outer2)])
    outer1, length1 = PerirhizalPython.to_range_(None, outer1, length1, 0., 2.)
    outer2, length2 = PerirhizalPython.to_range_(None, outer2, length2, 0., 2.)

    if i == 0:
        ax = axes[0, 0]
    else:
        ax = axes[0, 2]
    # ax.hist(outer_r, bins = 40, rwidth = 0.9)
    ax.hist([outer1, outer2], weights = [length1, length2], bins = 20, rwidth = 0.9, label = ["Top soil", "Sub soil"], stacked = True)
    ax.set_ylim(0, 440)
    ax.set_xlim(0., 2.)
    ax.set_title(titles[i])

    if i == 0:
        ax.set_ylabel("Root length [cm] ")

stats = np.array(stats)
for i, cell_number in enumerate(cell_numbers):
    print("\nSummary Spring Barley (min, max, median, mean, std):", cell_number)
    print("top soil", stats[i,:5])
    print("sub soil", stats[i, 5:])

axes[0, 1].axis('off')
ax.legend()

""" maize """
plant = "maize"
titles = []
titles.append("Maize 3D grid [1 cm]$^3$")
titles.append("Maize 2D grid [1 cm]$^2$")
titles.append("Maize 1D grid [1 cm]$^1$")
cell_numbers = []
cell_numbers.append([76, 16, 150])
cell_numbers.append([76, 1, 150])
cell_numbers.append([1, 1, 150])
dim_ = ["3D", "2D", "1D"]

stats = []
for i, dim in enumerate(dim_):

    outer1, outer2, seg2cell, length1, length2 = get_outer_radii(plant, type_str, dim)
    stats.append([np.min(outer1), np.max(outer1), np.median(outer1), np.mean(outer1), np.std(outer1),
                  np.min(outer2), np.max(outer2), np.median(outer2), np.mean(outer2), np.std(outer2)])
    outer1, length1 = PerirhizalPython.to_range_(None, outer1, length1, 0., 2.)
    outer2, length2 = PerirhizalPython.to_range_(None, outer2, length2, 0., 2.)

    ax = axes[1, i]
    # ax.hist(outer_r, bins = 40, rwidth = 0.9)
    ax.hist([outer1, outer2], weights = [length1, length2], bins = 20, rwidth = 0.9, label = ["Top soil", "Sub soil"], stacked = True)
    ax.set_ylim(0, 6250)
    # ax.set_xlim(0, 4)

    ax.set_title(titles[i])

    if i == 0:
        ax.set_ylabel("Root length [cm] ")
    ax.set_xlabel("Perirhizal outer radius [cm], 20 bins")

stats = np.array(stats)
for i, cell_number in enumerate(cell_numbers):
    print("\nSummary Maize (min, max, median, mean, std):", cell_number)
    print("top soil", stats[i,:5])
    print("sub soil", stats[i, 5:])

plt.tight_layout()
plt.savefig('outer_radii_grids_' + type_str + '.png')
plt.show()

print("fin")
