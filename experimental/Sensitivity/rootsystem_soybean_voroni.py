
""" 
Analyse and choose (modify) soybean system parameters 
"""
import sys; sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import plantbox as pb
import visualisation.vtk_plot as vp
import scenario_setup as scenario

from scipy.spatial import ConvexHull, Voronoi, voronoi_plot_2d
import vtk
from vtk import vtkXMLUnstructuredGridWriter, vtkUnstructuredGrid, vtkPoints, vtkDoubleArray

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

# typical domain for soybean
soil_, table_name, min_b, max_b, cell_number, area, Kc = scenario.soybean(0)
min_b = [-10., -10., -30.]  # data from INARI
max_b = [10., 10., 0.]
simtime = 5  # between 75-100 days

# Open plant and root parameter from a file
rs = pb.RootSystem()
path = "../../../CPlantBox/modelparameter/rootsystem/"
name = "Glycine_max_Moraes2020_opt2"
rs.readParameters(path + name + ".xml")

srp = rs.getOrganRandomParameter(pb.OrganTypes.seed)  # print(srp[0])
# srp[0].firstB = 1.e9  # 3
# srp[0].delayB = 3
src = srp[0].maxB

rrp = rs.getOrganRandomParameter(pb.OrganTypes.root)
for p in rrp:
    print("\nSubType", p.subType)
    print("Radius", p.a)
    print("dx", p.dx, "cm")
    print("theta", p.theta, "cm")
    print("lmax", p.lmax, "cm")
    print("changed to 0.5 cm to be faster...")
    p.dx = 0.5  # probably enough

p = np.array([1.* 2 ** x for x in np.linspace(-2., 2., 9)])

s, soil = scenario.create_soil_model(soil_, min_b, max_b, cell_number, type = 1)  # , times = x_, net_inf = y_
xml_name = "data/" + name + "_modified" + ".xml"  # root growth model parameter file
mods = {"lmax145":1., "lmax2":1., "lmax3":1., "theta45":1.5708, "r145":1., "r2":1., "a":1., "src":src}
r = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, xml_name, stochastic = False, mods = mods)  # pass parameter file for dynamic growth

kr = 1.e-4
kx = 1.e-3
scenario.init_conductivities_const(r, kr, kx)

rs = r.rs

# Initialize
print()
bounds = np.array(max_b) - (min_b)
rs.setGeometry(pb.SDF_PlantBox(*bounds))
rs.initialize()

# Simulate
rs.simulate(simtime, True)

vp.plot_roots(rs, "subType")
rs.write("rootsystem.vtp")
# Analyse
ana = pb.SegmentAnalyser(rs)
# mana.addAge(simtime)

orders = np.array(rs.getParameter("subType"))
print("\nnumber of roots", len(rs.getRoots()))
print("types", np.sum(orders == 1), np.sum(orders == 2), np.sum(orders == 3), np.sum(orders == 4), np.sum(orders == 5))
print("number of nodes", len(ana.nodes))
print("number of segments", len(ana.segments))
print("volume", np.sum(ana.getParameter("volume")), "cm3")
print("surface", np.sum(ana.getParameter("surface")), "cm2")
print("Krs", r.get_krs(simtime)[0], "cm2/day")

nodes = ana.nodes
nodes = np.array([[n.x, n.y, n.z] for n in nodes])
print(nodes.shape)


def voronoi_volumes(points):
    v = Voronoi(points)
    vol = np.zeros(v.npoints)
    for i, reg_num in enumerate(v.point_region):
        indices = v.regions[reg_num]
        if -1 in indices:  # some regions can be opened
            vol[i] = np.inf
        else:
            vol[i] = ConvexHull(v.vertices[indices]).volume
    return vol


points = nodes
vol = voronoi_volumes(points)
# print(vol)
vol0 = vol.copy()
vol0[vol == np.inf] = 0.
vol0[vol > 1.e3] = 0.
print("Total volume ", np.prod(bounds), "cm3")
print("Summed cells", np.sum(vol0[18:]), "cm3")

# Compute the Voronoi diagram
vor = Voronoi(points)

# Create a VTK unstructured grid
grid = vtkUnstructuredGrid()

# ana.mapPeriodic(w[0], w[1])
# ana.write("results/soybean_periodic.vtp")

# fig = plt.figure()
# ax = fig.add_subplot(111, projection = '3d')
#
# # # Plot the quader
# # for s, e in zip(quader, np.roll(quader, -1, axis = 0)):
# #     ax.plot([s[0], e[0]], [s[1], e[1]], [s[2], e[2]], 'k-', linewidth = 2)
#
# # Plot the Voronoi diagram
# for i, r in enumerate(vor.regions):
#     if -1 not in r:
#         polygon = [vor.vertices[j] for j in r]
#         ax.fill(*zip(*polygon), alpha = 0.3)
#
# # Set the axis limits
# ax.set_xlim([0, 1])
# ax.set_ylim([0, 1])
# ax.set_zlim([0, 1])
#
# plt.show()
#

ver = vor.vertices  # Veronoi vertices

# print("vertices")
# print(ver)
# print("point_region", vor.point_region)

# plt.plot(points[:, 0], points[:, 1], "b*", label = "points")
# plt.plot(ver[:, 0], ver[:, 1], "r*", label = "voronoi points")
# plt.legend()
# plt.show()

grid = vtkUnstructuredGrid()  # Create a VTK unstructured grid

points_array = vtkPoints()  # Add the Voronoi vertices as points
for v in ver:
    for j in range(0, 3):
        v[j] = min(v[j], max_b[j])
        v[j] = max(v[j], min_b[j])
    points_array.InsertNextPoint(v)
grid.SetPoints(points_array)
print("voronoi nodes: ", ver.shape)

# print("vor.regions", vor.regions)
cc = 0
for region in vor.regions:
    if len(region) > 0 and (-1 not in region):
        cc += 1
        id_array = vtk.vtkIdList()
        id_array.InsertNextId(len(region))
        # print("added an id array with length", len(region))
        for vertex_index in region:
            id_array.InsertNextId(vertex_index)
        grid.InsertNextCell(vtk.VTK_CONVEX_POINT_SET, id_array)

cell_id = vtk.vtkDoubleArray()
cell_id.SetName("cell_id")
cell_id.SetNumberOfValues(cc)
vc = 0
for j in range(0, cc):
    cell_id.SetValue(vc, vc)
    vc += 1
celldata = grid.GetCellData()
celldata.AddArray(cell_id)

cell_id = vtk.vtkDoubleArray()
cell_id.SetName("volumes")
cell_id.SetNumberOfValues(cc)
vc = 0
for j, region in enumerate(vor.regions):
    if len(region) > 0 and (-1 not in region):
        cell_id.SetValue(vc, vol0[j])
        vc += 1
celldata = grid.GetCellData()
celldata.AddArray(cell_id)

print("alive")

# Write the VTK unstructured grid to a ParaView VTU file
writer = vtkXMLUnstructuredGridWriter()
writer.SetFileName('rootsystem.vtu')
writer.SetInputData(grid)
writer.Write()
print("fin")

# quader_polygons = []
# for i, region in enumerate(vor.regions):
#     if -1 in region:  # Skip the infinite region
#         continue
#     polygon = [vor.vertices[j] for j in region]
#     intersection = False
#     for s, e in zip(quader, np.roll(quader, -1, axis = 0)):
#         # Check if the polygon intersects with any of the quader faces
#         if VoronoiEdgeIntersection(polygon, s, e):
#             intersection = True
#             break
#     if intersection:
#         quader_polygons.append(polygon)
#
#
# def VoronoiEdgeIntersection(polygon, s, e):
#     for p1, p2 in zip(polygon, np.roll(polygon, -1, axis = 0)):
#         # Check if the polygon edge intersects with the quader edge
#         if LineSegmentIntersection(p1, p2, s, e):
#             return True
#     return False
#
#
# def LineSegmentIntersection(p1, p2, p3, p4):
#     # Compute the direction vectors of the line segments
#     d1 = p2 - p1
#     d2 = p4 - p3
#     # Compute the cross products of the direction vectors
#     c1 = np.cross(d1, p1 - p3)
#     c2 = np.cross(d1, d2)
#     # Check if the line segments are parallel
#     if np.allclose(c2, 0):
#         return False
#     # Compute the intersection point
#     t = np.cross(p3 - p1, d2) / np.dot(d1, d2)
#     intersection = p1 + t * d1
#     # Check if the intersection point lies within both line segments
#     return (0 <= t <= 1) and (0 <= np.dot(intersection - p3, d2) <= np.dot(d2, d2))

