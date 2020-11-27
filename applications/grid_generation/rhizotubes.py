import distmesh as dm
import numpy as np
import sdf

from vtk_tools import *
import vtk

by = 1.
bx = .5
bz = 1.5


def rhizo_tubes():
    """ Creates the rhizo tubes geometry
    """
    r = 3.2 / 100.  # tube radius
    l = bx  # cm tube length
    rhizotube = sdf.Plant_Container(r, r, l, False)
    rhizoX = sdf.Rotate_Translate(rhizotube, np.array([l, 0., 0.]), 90., 1)
    y_ = np.array([ 35, 45, 55, 65, 75, 85 ]) / 100
    z_ = np.array([ -10, -20, -40, -60, -80, -120 ]) / 100.
    rhizotubes_ = []
    for i in range(0, len(y_)):
        newtube = sdf.Rotate_Translate(rhizoX, np.array([0, y_[i], z_[i]]), 0., 0)
        rhizotubes_.append(newtube)
    return sdf.Union(rhizotubes_)


bbox = np.array([0, 0, -bz, bx, by, 0])
container = sdf.Box(bbox)
tubes = rhizo_tubes()
geom = sdf.Difference(container, tubes)
fh = sdf.Edge_Length([tubes], 0.2, 1, 0.1)
p, t = dm.distmeshnd(geom.f, fh.f, 0.005, bbox)

print("rebuild")
p, t = rebuild_grid(p, t)
print("snap")
p = snap_to_box(p, bbox, 1e-5)
name = "rhizo.msh"

p[:, 2] += 2  # shift from 0 to  +200

points = vtkPoints(p)  # make the unstructured grid
cells = vtkCells(t)
grid = vtk.vtkUnstructuredGrid()
grid.SetPoints(points)
grid.SetCells(vtk.VTK_TETRA, cells)
celldata = np.hstack((np.ones((grid.GetNumberOfCells(), 1)), 10 * np.ones((grid.GetNumberOfCells(), 1))))

# write vtu
writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetFileName(name[0:-4] + ".vtu");
writer.SetInputData(grid);
writer.Write();
write_msh(name, grid, celldata)

print()
print("Points: ", p.shape)
print("Triangles", t.shape)
print()
grid_quality(p, t)

print("done")
