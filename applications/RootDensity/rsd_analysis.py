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

""" parameters """
simtime = 87.5  # between 75-100 days

""" setup """
soil_, table_name, min_b, max_b, cell_number, area, Kc = scenario.soybean(0)  # 0 = envirotype
xml_name = "data/Glycine_max_Moraes2020_opt2_modified.xml"  # root growth model parameter file
cell_number = [76, 4, 100]
# cell_number = [38, 2, 50]
# cell_number = [19, 1, 25]

# soil_, table_name, min_b, max_b, cell_number, area, Kc = scenario.maize(0)  # 0 = envirotype
# xml_name = "data/Zeamays_synMRI_modified.xml"  # root growth model parameter file
# cell_number = [76, 16, 200]
# # cell_number = [38, 8, 100]
# # cell_number = [19, 4, 50]

width = max_b - min_b

s, soil = scenario.create_soil_model(soil_, min_b, max_b, cell_number, type = 1)
# sra_table_lookup = sra.open_sra_lookup("data/" + table_name)

r_ = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, xml_name)  # pass parameter file for dynamic growth
# scenario.init_conductivities_const_growth(r)
# scenario.init_lupine_conductivities(r)

r = r_.rs  # throw away (TODO ahve to change setup anyway...)

r.initialize(False)
r.simulate(simtime, True)

sn = np.prod(cell_number)
peri = PerirhizalPython(r)
sd = peri.get_density("surface")

grid = vp.uniform_grid(min_b, max_b, cell_number)
cell_sd = vtk.vtkDoubleArray()
cell_sd.SetName("surface_density")
cell_sd.SetNumberOfValues(sn)
for j in range(0, sn):
    cell_sd.SetValue(j, min(sd[j], 1.))
celldata = grid.GetCellData()
celldata.AddArray(cell_sd)

outer_radii = peri.get_outer_radii_voronoi()
print("outer_radii", np.min(outer_radii), np.max(outer_radii), "median", np.median(outer_radii), "mean", np.mean(outer_radii), np.std(outer_radii))

ana = pb.SegmentAnalyser(r.mappedSegments())
print(len(ana.segments))
print(len(outer_radii))
# outer_radii = np.minimum(outer_radii, 1.)  # limit for visualisation
ana.addData("outer_r", outer_radii)

# vp.plot_mesh(grid, "surface_density")
# vp.plot_mesh_cuts(grid, "surface_density")
# vp.plot_roots(ana, "outer_r")
vp.plot_roots_and_mesh(ana, "outer_r", grid, "surface_density", True, width[0], width[1])

# plt.hist(outer_radii, bins = 100, rwidth = 0.9)
# plt.show()

print("fin")
