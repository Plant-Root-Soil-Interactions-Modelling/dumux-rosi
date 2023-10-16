""" 
    plots histograms of the outer radius of the perirhizal zone for three soil horizons using Voronoi diagrams
    
    copy & paste of plot_outerradii_horizons.py, with other caption and L57 changed
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


def get_outer_radii_horizons(rootsystem):

    """ setup """
    if rootsystem == "Soybean":
        soil_, table_name, min_b, max_b, _, area, Kc = scenario.soybean(0)  # 0 = envirotype
        xml_name = "data/Glycine_max_Moraes2020_opt2_modified.xml"  # root growth model parameter file
        cell_number = [38, 2, 150]  # (2cm)^3
        cell_number = [1, 1, 150]  # (2cm)^3
#        cell_number = [1, 1, 1]  # (2cm)^3
        simtime = 42  # 87.5
    elif rootsystem == "Maize":
        soil_, table_name, min_b, max_b, _, area, Kc = scenario.maize(0)  # 0 = envirotype
        xml_name = "data/Zeamays_synMRI_modified.xml"  # root growth model parameter file
        simtime = 56  #
        # cell_number = [38, 8, 159]  # (2cm)^3
        # cell_number = [1, 1, 159]  # (2cm)^3
        cell_number = [1, 1, 1]  # (2cm)^3
    elif rootsystem == "Spring Barley":
        soil_, table_name, min_b, max_b, _, area, Kc = scenario.springbarley(0)  # 0 = envirotype
        xml_name = "data/spring_barley_CF12.xml"  # root growth model parameter file
        simtime = 49  # 95
        # cell_number = [7, 2, 75]  # (2cm)^3
        cell_number = [1, 1, 1]  # (2cm)^3
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
    # outer_radii = peri.get_outer_radii_voronoi()
    outer_radii = peri.get_outer_radii_bounded_voronoi()  ##############################################################################################

    print()
    print("r.nodes", len(r.nodes))
    print("outer_radii", len(outer_radii))
    print(np.argwhere(np.isnan(outer_radii)))
    print()

    print()
    print("outer_radii", np.nanmin(outer_radii), np.nanmax(outer_radii), "median", np.nanmedian(outer_radii), "mean", np.nanmean(outer_radii), np.nanstd(outer_radii))
    print()

    print("outer_radii", outer_radii.shape)
    # for r in outer_radii:
    #     print(r)

    ana = pb.SegmentAnalyser(r.mappedSegments())
    ana.addData("outer_r", outer_radii)
    ana.addData("length", ana.getParameter("length"))  # to make sure that we handle it exactly the same way

    # organic_layer = pb.SDF_Cuboid(pb.Vector3d([-1e6, -1e6, -organic]), pb.Vector3d([1e6, 1e6, max_b[2]]))
    topsoil_layer = pb.SDF_Cuboid(pb.Vector3d([-1e6, -1e6, -topsoil]), pb.Vector3d([1e6, 1e6, max_b[2]]))
    subsoil_layer = pb.SDF_Cuboid(pb.Vector3d([-1e6, -1e6, min_b[2]]), pb.Vector3d([1e6, 1e6, -topsoil]))
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
    print(len(outer_radii), len(outer1) + len(outer2))

    return outer1, outer2, length1, length2


""" parameters """
# see https://www.soils4teachers.org/soil-horizons/
# organic = 6  # 2 * 2.54
topsoil = 30  # 10 * 2.54
subsoil = 150  # 30 * 2.54

rootsystem = "Spring Barley"  # Maize, Soybean, Spring Barley

# fig, axes = plt.subplots(3, 1, figsize = (10, 18))
fig, axes = plt.subplots(1, 1, figsize = (9, 8))
outer1, outer2, length1, length2 = get_outer_radii_horizons(rootsystem)

outer1, length1 = PerirhizalPython.to_range_(None, outer1, length1, 0., 2.)
outer2, length2 = PerirhizalPython.to_range_(None, outer2, length2, 0., 2.)

axes.hist([outer1, outer2], weights = [length1, length2], bins = 40, rwidth = 0.9, label = ["topsoil", "subsoil"], stacked = True)  # outer0,
#  axes.hist([outer1, outer2], bins = 40, rwidth = 0.9, label = ["topsoil", "subsoil"], stacked = True)  # outer0,

# axes.set_ylim([0., 4000.])
axes.legend()
axes.set_xlabel("Perirhizal outer radius [cm] using Voronoi diagrams, 40 bins")
axes.set_ylabel("Root length [cm] ")
axes.set_title(rootsystem)
plt.tight_layout()
plt.savefig('hist_vor_' + rootsystem + '.png')

plt.show()

print("fin")
