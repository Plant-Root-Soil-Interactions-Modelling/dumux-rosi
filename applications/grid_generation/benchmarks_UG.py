import sys; sys.path.append("../../../python/modules/"); sys.path.append("../../../../CPlantBox/"); sys.path.append("../../../CPlantBox/src/python_modules")

from distmeshnd import *
import sdf
from vtk_tools import *
import numpy as np


def huniform(p):
    """Implements the trivial uniform mesh size function h=1."""
    return np.ones(p.shape[0])


def fh_soil_layer(p, layerZ, hmin = 0.5, hmax = 2., slope = 1):
    """ Higher resolution around one or more layers """
    for l in layerZ:
        h1 = np.maximum(slope * (np.abs(p[:, 2] - l)), hmin)
    return np.minimum(h1, hmax)


def tet_mid(tt, p):
    """ returns the midpoints of the tetraedras """
    midp = np.zeros((tt.shape[0], 3))
    for i, t in enumerate(tt):
        midp[i, :] = 0.25 * (p[t[0], :] + p[t[1], :] + p[t[2], :] + p[t[3], :])
    return midp


celldata = np.zeros((0, 1))

#
# BENCHMARK 1
#
# 9*9*199 = 16119
# r^2*pi*199 = 12659.5 ( r = 4.5 cm )

box = (0, 0, 0, 0.1, 0.2, 2)  # b1
fBox = sdf.Box(box)
fun = lambda p : fBox.f(p)
cbox = (-0.1, -0.1, 0, 0.1, 0.1, 2)
cylinder = sdf.Plant_Container(0.1, 0.1, 2)

h0 = 0.029  # intial edge lenght (m), (controls dof)
p, t = distmeshnd(cylinder.f, huniform, h0, np.array(cbox), None, 500)
p, t = rebuild_grid(p, t)
p = snap_to_box(p, cbox, 1e-5)
name = "b1_ug.msh"  # 9k dof

# fh1 = lambda p: fh_soil_layer(p, [2,1.5], 0.4, 1, 3)
# p, t = distmeshnd(cylinder.f, fh1, 0.01, np.array(cbox), pfix = np.array([ [0.1,0,0],[-0.1,0,0],[0,0.1,0],[0,-0.1,0], [0.1,0,2],[-0.1,0,2],[0,0.1,2],[0,-0.1,2] ]))
# name = "b1_ug2.msh" # 6K dof

celldata = np.zeros((t.shape[0], 1))
tetmid = tet_mid(t, p)
for i, t_ in enumerate(tetmid):
    z = t_[2]
    if z >= 1.5:
        celldata[i, 0] = 1
    else:
        celldata[i, 0] = 2

# #
# # BENCHMARK 2
# #
# cbox2 = (-0.1, -0.1, 0, 0.1, 0.1, .54)
# cylinder2 = sdf.Plant_Container(0.1, 0.1, .54)
#
# p, t = distmeshnd(cylinder2.f, dm.huniform, 0.01, np.array(cbox2))
# p, t = rebuild_grid(p, t)
# p = snap_to_box(p, cbox2, 1e-5)
# name = "b2_ug.msh"  # 17k
#
# # fh2 = lambda p: fh_soil_layer(p, [0.54], 0.3, 1, 3)
# # p, t = distmeshnd(cylinder2.f, fh2, 0.0075, np.array(cbox2), pfix = np.array([ [0.1,0,0],[-0.1,0,0],[0,0.1,0],[0,-0.1,0], [0.1,0,.54],[-0.1,0,.54],[0,0.1,.54],[0,-0.1,.54] ]))
# # p, t = rebuild_grid(p,t)
# # p = snap_to_box(p, cbox2, 1e-5)
# # name = "b2_ug2.msh" # 11K dof

#
# Make the unstructured grid
#
points = vtk_points(p)
cells = vtk_cells(t)
grid = vtk.vtkUnstructuredGrid()
grid.SetPoints(points)
grid.SetCells(vtk.VTK_TETRA, cells)
if celldata.shape[0] > 0:
    grid.GetCellData().SetScalars(vtk_data(celldata))

#
print()
print("Points: ", p.shape)
print("Triangles", t.shape)
print()

# write vtu
writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetFileName(name[0:-4] + ".vtu");
writer.SetInputData(grid);
writer.Write();

write_msh(name, grid)

# vtkInterface.Plot(grid)

print("done")
