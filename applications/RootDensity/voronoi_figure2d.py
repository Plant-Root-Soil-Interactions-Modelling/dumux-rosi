""" figures for the icons (2D)"""
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

r = r_.rs  # throw away (TODO ahve to change setup anyway...)

r.initialize(False)
r.simulate(simtime, False)

r.write("rootsystem_geometry_large.py")
ana = pb.SegmentAnalyser(r.mappedSegments())
ana.mapPeriodic(width[0], width[1])
ana.write("rootsystem_large.vtp")

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

outer_radii = peri.get_outer_radii_voronoi()
# outer_radii = peri.get_outer_radii("surface")

print()
print("number of outer", outer_radii.shape)
print("open regions", np.sum(outer_radii == -1), np.sum(outer_radii == -1) / outer_radii.shape[0])
print("cell has point outside domain", np.sum(outer_radii == 0), np.sum(outer_radii == 0) / outer_radii.shape[0])
print()

print("outer_radii", np.min(outer_radii), np.max(outer_radii), "median", np.median(outer_radii), "mean", np.mean(outer_radii), np.std(outer_radii))
# vp.plot_mesh(grid, "surface_density")
# vp.plot_mesh_cuts(grid, "surface_density")
# vp.plot_roots(ana, "outer_r")

""" outer_r 3d visualisation """
ana = pb.SegmentAnalyser(r.mappedSegments())
outer_radii = -np.minimum(outer_radii, 1.)  # limit for visualisation
ana.addData("outer_r", outer_radii)
vp.plot_roots_and_mesh(ana, "subType", grid, "surface_density", True, width[0], width[1])  # outer_r

""" nice histogram """
fig, axes = plt.subplots(1, 1, figsize = (8, 8))
rr = peri.to_range_(outer_radii, 0., 2.)
axes.hist(rr, bins = 40, rwidth = 0.9)
plt.show()

# Write the VTK unstructured grid to a ParaView VTU file
domain = pb.SDF_Cuboid(pb.Vector3d(min_b), pb.Vector3d(max_b))
grid = peri.get_voronoi_mesh(domain)
writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetFileName('rootsystem_mesh_large.vtu')
writer.SetInputData(grid)
writer.Write()

print("fin")
