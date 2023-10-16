import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull, Voronoi, voronoi_plot_2d
import vtk
from vtk import vtkXMLUnstructuredGridWriter, vtkUnstructuredGrid, vtkPoints, vtkDoubleArray, vtkIdList

rng = np.random.default_rng()
points_ = rng.random((30, 3))  # 30 random points in 2-D
points = []
for p in points_:
    if (p > -1).all() and (p < 1).all():
        points.append(p)

points = np.array(points)
print("points", points.shape)

vor = Voronoi(points)
ver = vor.vertices  # Veronoi vertices

# print("vertices")
# print(ver)
# print("point_region", vor.point_region)

plt.plot(points[:, 0], points[:, 1], "b*", label = "points")
plt.plot(ver[:, 0], ver[:, 1], "r*", label = "voronoi points")
plt.legend()
plt.show()

grid = vtkUnstructuredGrid()  # Create a VTK unstructured grid

points_array = vtkPoints()  # Add the Voronoi vertices as points
for v in ver:
    points_array.InsertNextPoint(v)
grid.SetPoints(points_array)
print("voronoi nodes: ", ver.shape)

# print("vor.regions", vor.regions)
cc = 0
for region in vor.regions:
    if len(region) > 0 and (-1 not in region):
        cc += 1
        id_array = vtkIdList()
        id_array.InsertNextId(len(region))
        print("added an id array with length", len(region))
        for vertex_index in region:
            id_array.InsertNextId(vertex_index)
        grid.InsertNextCell(vtk.VTK_CONVEX_POINT_SET, id_array)

# convex = []
# for i, reg_num in enumerate(vor.point_region):
#     indices = vor.regions[reg_num]
#     if -1 in indices:  # some regions can be opened
#         convex.append(None)
#     else:
#         convex.append(ConvexHull(ver[indices]))

# cc = 0
# for j, c in enumerate(convex):
#     if c:
#         cc += 1
#         faceId = vtkIdList()
#         # faceId.InsertNextId(len(c.simplices))  # number of cells
#         print("convex shape", len(c.simplices), "simplices", convex[j].volume)
#         for simplex in c.simplices:
#             print(len(simplex))
#             faceId.InsertNextId(3)
#             [faceId.InsertNextId(i) for i in simplex]    cell_id.SetValue(vc, vc)
    vc += 1
#         grid.InsertNextCell(vtk.VTK_POLYHEDRON, faceId)
#         break
#         # for vertices_id in c.simplices:
#         #     print(vertices_id, ", ", end = "")
#         # print()

cell_id = vtk.vtkDoubleArray()
cell_id.SetName("cell_id")
cell_id.SetNumberOfValues(cc)
vc = 0
for j in range(0, cc):
    cell_id.SetValue(vc, vc)
    vc += 1
celldata = grid.GetCellData()
celldata.AddArray(cell_id)

print("alive")

# Write the VTK unstructured grid to a ParaView VTU file
writer = vtkXMLUnstructuredGridWriter()
writer.SetFileName('voronoi.vtu')
writer.SetInputData(grid)
writer.Write()
print("fin")
