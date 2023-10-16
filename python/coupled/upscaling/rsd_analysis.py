""" 
    3d surface densities
"""
import sys; sys.path.append("../../modules"); sys.path.append("../../../build-cmake/cpp/python_binding/");
sys.path.append("../../../../CPlantBox");  sys.path.append("../../../../CPlantBox/src")

import plantbox as pb
from functional.Perirhizal import *
import visualisation.vtk_plot as vp
import scenario_setup as scenario

import vtk
import numpy as np
import matplotlib.pyplot as plt

""" parameters """

dim = "3D"  # 1D, 3D
plant = "springbarley"  # soybean, maize, springbarley
# min_b, max_b, cell_number = scenario.soybean_(dim)
min_b, max_b, cell_number = scenario.springbarley_(dim)
soil = "hydrus_loam"  #  hydrus_loam, hydrus_clay, hydrus_sand
outer_method = "surface"  # voronoi, length, surface, volume
initial = -200  # cm total potential

r_, rho_, rs_age, trans, wilting_point, soil, s, sra_table_lookup, mapping = scenario.set_scenario(plant, dim, initial, soil, outer_method)
r = r_.rs  # throw away (TODO ahve to change setup anyway...)

segs = r.segments
nodes = r.nodes

print("\nplant")
print("rs_age", rs_age)
print()
params = r.getRootRandomParameter()
for p in params:
    print("SubType", p.subType)
    print("radius", p.a)
    print("growth rate", p.r)
    print()

print()

sn = np.prod(cell_number)
peri = PerirhizalPython(r)

""" add 3d soil surface density """
sd = peri.get_density("surface")
sd = -np.minimum(sd, 1.1)  # limit for visualisation
grid = vp.uniform_grid(min_b, max_b, cell_number)
cell_sd = vtk.vtkDoubleArray()
cell_sd.SetName("surface_density")
cell_sd.SetNumberOfValues(sn)
for j in range(0, sn):
    cell_sd.SetValue(j, sd[j])
celldata = grid.GetCellData()
celldata.AddArray(cell_sd)

# outer_radii = peri.get_outer_radii_voronoi()
outer_radii = peri.get_outer_radii("surface")

print("outer_radii", np.min(outer_radii), np.max(outer_radii), "median", np.median(outer_radii), "mean", np.mean(outer_radii), np.std(outer_radii))

ana = pb.SegmentAnalyser(r.mappedSegments())
# outer_radii = -np.minimum(outer_radii, 1.)  # limit for visualisation
ana.addData("outer_r", outer_radii)
ana.addHydraulicConductivities(r_.params, rs_age, 0.1, 0.014)
ana.addAge(rs_age)
ana.addCellIds(r)

# vp.plot_mesh(grid, "surface_density")
vp.plot_mesh_cuts(grid, "surface_density")
print()
print("age", rs_age, "nodes", len(nodes))
print()

# vp.plot_roots(ana, "age")
# vp.plot_roots(ana, "kr")
# vp.plot_roots(ana, "kx")
vp.plot_roots(ana, "subType")

ana.write("results/" + plant + ".vtp")

# width = max_b - min_b
# vp.plot_roots_and_mesh(ana, "outer_r", grid, "surface_density", True, width[0], width[1])

# plt.hist(outer_radii, bins = 40, rwidth = 0.9, align = 'mid')
# plt.show()

print("fin")
