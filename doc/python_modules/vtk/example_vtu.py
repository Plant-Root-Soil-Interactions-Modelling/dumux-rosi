import sys; sys.path.append("../../../python/modules/"); sys.path.append("../../../../CPlantBox/")

import vtk_plot as vp
import vtk_tools as vt

name = "soybean_Honly-00001"

pd = vp.read_vtu(name + ".vtu")  # or read_vtu

# np_data(polydata, data_index = 0, cell = None):
vp.plot_roots_and_soil_files(name, "pressure head")

print("fin")
