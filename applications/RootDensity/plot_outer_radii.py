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

""" parameters """

type_str = "surface"

# see https://www.soils4teachers.org/soil-horizons/
organic = 6  # 2 * 2.54
topsoil = 26  # 10 * 2.54
subsoil = 76  # 30 * 2.54

soil_, table_name, min_b, max_b, cell_number, area, Kc = scenario.soybean(0)  # 0 = envirotype
simtime = 87.5  # between 75-100 days
xml_name = "data/Glycine_max_Moraes2020_opt2_modified.xml"  # root growth model parameter file
cell_number = [76, 4, 100]
# cell_number = [38, 2, 50]
# cell_number = [19, 1, 25]

# soil_, table_name, min_b, max_b, cell_number, area, Kc = scenario.maize(0)  # 0 = envirotype
# xml_name = "data/Zeamays_synMRI_modified.xml"  # root growth model parameter file
# simtime = 95.5  # between 75-100 days
# # cell_number = [76, 16, 200]
# cell_number = [38, 8, 100]  # (2cm)^3
# # cell_number = [19, 4, 50]

width = max_b - min_b

s, soil = scenario.create_soil_model(soil_, min_b, max_b, cell_number, type = 1)
r_ = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, xml_name)  # pass parameter file for dynamic growth
r = r_.rs
r.initialize(True)
r.simulate(simtime, True)

""" analysis """
peri = PerirhizalPython(r)
outer_radii = peri.get_outer_radii(type_str)
outer_radii = np.minimum(outer_radii, 2.)
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

# vp.plot_mesh(grid, "surface_density")
# vp.plot_mesh_cuts(grid, "surface_density")
# vp.plot_roots(ana, "outer_r")
# vp.plot_roots_and_mesh(ana, "outer_r", grid, "surface_density", True, width[0], width[1])
outer0 = ana0.data["outer_r"]
outer1 = ana1.data["outer_r"]
outer2 = ana2.data["outer_r"]

print(len(outer_radii), len(outer0) + len(outer1) + len(outer2))

fig, axes = plt.subplots(2, 1, figsize = (18, 10))

axes[0].hist(outer_radii, bins = 40, rwidth = 0.9)
axes[1].hist([outer0, outer1, outer2], bins = 40, rwidth = 0.9, label = ["organic", "topsoil", "subsoil"], stacked = True)
axes[1].legend()
plt.show()

print("fin")
