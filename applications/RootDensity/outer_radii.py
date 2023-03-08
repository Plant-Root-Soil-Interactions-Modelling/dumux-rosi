"""small example"""
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


def get_outer_radii(rootsystem, type, cell_number):

    """ setup """
    if rootsystem == "soybean":
        soil_, table_name, min_b, max_b, _, area, Kc = scenario.soybean(0)  # 0 = envirotype
        xml_name = "data/Glycine_max_Moraes2020_opt2_modified.xml"  # root growth model parameter file
        simtime = 87.5
    elif rootsystem == "maize":
        soil_, table_name, min_b, max_b, _, area, Kc = scenario.maize(0)  # 0 = envirotype
        xml_name = "data/Zeamays_synMRI_modified.xml"  # root growth model parameter file
        simtime = 95

    width = max_b - min_b
    s, soil = scenario.create_soil_model(soil_, min_b, max_b, cell_number, type = 1)
    r_ = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, xml_name)  # pass parameter file for dynamic growth

    """ simulation and post processing """
    r = r_.rs  # rename
    r.initialize(False)
    r.simulate(simtime, False)
    sn = np.prod(cell_number)
    peri = PerirhizalPython(r)
    outer_radii = peri.get_outer_radii(type)
    # sd = peri.get_density("surface")

    return outer_radii


fig, axes = plt.subplots(3, 2, figsize = (10, 18))

""" maize """
cell_numbers = []
cell_numbers.append([76, 16, 200])
cell_numbers.append([38, 8, 100])
cell_numbers.append([19, 4, 50])
cell_numbers.append([1, 1, 200])
cell_numbers.append([1, 1, 100])
cell_numbers.append([1, 1, 50])

stats = [["min", "max", "median", "mean", "std"]]
for i, cell_number in enumerate(cell_numbers):
    outer_r = get_outer_radii("maize", "surface", cell_number)
    outer_r = np.minimum(outer_r, 2)
    ax = axes[i % 3, i // 3]
    ax.hist(outer_r, bins = 40, rwidth = 0.9)
    # ax.set_xlim(0, 4)
    ax.set_ylim(0, 15000)
    ax.set_title("resolution " + str(cell_number))
    stats.append([np.min(outer_r), np.max(outer_r), np.median(outer_r), np.mean(outer_r), np.std(outer_r)])

print(stats[0])
for i, cell_number in enumerate(cell_numbers):
    print("Summary", cell_number, ":")
    print(stats[i + 1])

plt.title("Maize")
plt.tight_layout()
plt.show()

# """ soybean """
# cell_numbers = []
# cell_numbers.append([76, 4, 200])
# cell_numbers.append([38, 2, 100])
# cell_numbers.append([19, 1, 50])
# cell_numbers.append([1, 1, 200])
# cell_numbers.append([1, 1, 100])
# cell_numbers.append([1, 1, 50])
#
# stats = [["min", "max", "median", "mean", "std"]]
# for i, cell_number in enumerate(cell_numbers):
#     outer_r = get_outer_radii("soybean", "surface", cell_number)
#     outer_r = np.minimum(outer_r, 2)
#     ax = axes[i % 3, i // 3]
#     ax.hist(outer_r, bins = 40, rwidth = 0.9)
#     # ax.set_xlim(0, 4)
#     ax.set_ylim(0, 3000)
#     ax.set_title("resolution " + str(cell_number))
#     stats.append([np.min(outer_r), np.max(outer_r), np.median(outer_r), np.mean(outer_r), np.std(outer_r)])
#
# print(stats[0])
# for i, cell_number in enumerate(cell_numbers):
#     print("Summary", cell_number, ":")
#     print(stats[i + 1])
#
# plt.title("Soybean")
# plt.tight_layout()
# plt.show()

print("fin")
