""" 
    plots histograms of the outer radius of the perirhizal zone for three soil horizons
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


def get_outer_radii_horizons(rootsystem, type_str):

    """ setup """
    if rootsystem == "Soybean":
        soil_, table_name, min_b, max_b, _, area, Kc = scenario.soybean(0)  # 0 = envirotype
        xml_name = "data/Glycine_max_Moraes2020_opt2_modified.xml"  # root growth model parameter file
        cell_number = [38, 2, 100]  # (2cm)^3
        simtime = 87.5
    elif rootsystem == "Maize":
        soil_, table_name, min_b, max_b, _, area, Kc = scenario.maize(0)  # 0 = envirotype
        xml_name = "data/Zeamays_synMRI_modified.xml"  # root growth model parameter file
        simtime = 95
        cell_number = [38, 8, 100]  # (2cm)^3
    else:
        print("get_outer_radii_horizons() unknown rootsystem name", rootsystem)
        raise

    """ simulate """
    s, soil = scenario.create_soil_model(soil_, min_b, max_b, cell_number, type = 1)
    r_ = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, xml_name)  # pass parameter file for dynamic growth
    r = r_.rs
    r.initialize(False)
    r.simulate(simtime, False)

    """ analysis """
    peri = PerirhizalPython(r)
    outer_radii = peri.get_outer_radii(type_str)
    print()
    print("outer_radii", np.min(outer_radii), np.max(outer_radii), "median", np.median(outer_radii), "mean", np.mean(outer_radii), np.std(outer_radii))
    print()
    ana = pb.SegmentAnalyser(r.mappedSegments())
    ana.addData("outer_r", outer_radii)
    organic_layer = pb.SDF_Cuboid(pb.Vector3d([-1e6, -1e6, -organic]), pb.Vector3d([1e6, 1e6, max_b[2]]))
    topsoil_layer = pb.SDF_Cuboid(pb.Vector3d([-1e6, -1e6, -topsoil]), pb.Vector3d([1e6, 1e6, -organic]))
    subsoil_layer = pb.SDF_Cuboid(pb.Vector3d([-1e6, -1e6, min_b[2]]), pb.Vector3d([1e6, 1e6, -topsoil]))
    ana0 = pb.SegmentAnalyser(ana)
    ana0.crop(organic_layer)
    ana0.pack()
    ana1 = pb.SegmentAnalyser(ana)
    ana1.crop(topsoil_layer)
    ana1.pack()
    ana2 = pb.SegmentAnalyser(ana)
    ana2.crop(subsoil_layer)
    ana2.pack()
    outer0 = ana0.data["outer_r"]
    outer1 = ana1.data["outer_r"]
    outer2 = ana2.data["outer_r"]
    print(len(outer_radii), len(outer0) + len(outer1) + len(outer2))
    return outer0, outer1, outer2


""" parameters """
# see https://www.soils4teachers.org/soil-horizons/
organic = 6  # 2 * 2.54
topsoil = 26  # 10 * 2.54
subsoil = 76  # 30 * 2.54 ############################################################ CHECK !!!!!

rootsystem = "Soybean"  # Maize, Soybean

# fig, axes = plt.subplots(3, 1, figsize = (10, 18))
fig, axes = plt.subplots(1, 1, figsize = (9, 8))
axes = [axes]

# types = ["length", "surface", "volume"]
types = ["length"]
for i in range(0, len(types)):
    outer0, outer1, outer2 = get_outer_radii_horizons(rootsystem, types[i])
    outer0 = PerirhizalPython.to_range_(None, outer0, 0., 2.)
    outer1 = PerirhizalPython.to_range_(None, outer1, 0., 2.)
    outer2 = PerirhizalPython.to_range_(None, outer2, 0., 2.)
    axes[i].hist([outer0, outer1, outer2], bins = 40, rwidth = 0.9, label = ["organic", "topsoil", "subsoil"], stacked = True)
    axes[i].legend()
    axes[i].set_xlabel("Perirhizal outer radius [cm] - proportional to " + types[i])
    axes[i].set_title(rootsystem)

plt.tight_layout()
plt.show()

print("fin")