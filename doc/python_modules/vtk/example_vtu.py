import sys; sys.path.append("../../../python/modules/"); sys.path.append("../../../../CPlantBox/")

import vtk_plot as vp
import vtk_tools as vt

"""
plot DuMux .vtu output

pull out data

make own PolyData grid, put in data, plot again
"""

name = "soybean_Honly-00001"

# Openn .vtu
pd = vp.read_vtu(name + ".vtu")
print(pd.GetBounds()) # xmin, xmax, ymin, ymax, zmin, zmax
print("Number of cells", vt.np_cells(pd).shape[0])

# Convert m -> cm
np_points = vt.np_points(pd)
points = vt.vtk_points(np_points*100) # m -> cm
pd.SetPoints(points)

# vp.plot_mesh(pd, "pressure head") # useless plot
vp.plot_mesh_cuts(pd, "pressure head", 5)

# np_data(polydata, data_index = 0, cell = None):

# vp.plot_roots_and_soil_files(name, "pressure head")

min_ = np.array([-18.5, -3, -150])
max_ = np.array([18.5, 3, 0.])
res_ = np.array([37, 6, 150])
rs.setRectangularGrid(pb.Vector3d(min_), pb.Vector3d(max_), pb.Vector3d(res_), True)  # cut and map segments

data = vt.np_data(pd,"pressure head", True) # we need cell data

# cplantbox 3d grid





print("fin")
