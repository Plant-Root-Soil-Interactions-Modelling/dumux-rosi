import numpy as np

import vtk_plot as vp
import vtk_tools as vt

""" Parameters """
min_b = [-4., -4., -15.]
max_b = [4., 4., 0.]
cell_number = [8, 8, 15]  # [32, 32, 60]  # [8, 8, 15]  # 32, 32, 60

soil_grid = vp.uniform_grid(np.array(min_b), np.array(max_b), np.array(cell_number))

soil_water_content = vt.vtk_data(list(range(0, np.prod(cell_number)))) # ordering... z, y, x 
soil_water_content.SetName("water content")
soil_grid.GetCellData().AddArray(soil_water_content)

soil_pressure = vt.vtk_data(np.ones((np.prod(cell_number),)))
soil_pressure.SetName("pressure head") # in macroscopic soil
soil_grid.GetCellData().AddArray(soil_pressure)

name = "water content"
soil_grid.GetCellData().SetActiveScalars(name)


meshActor, meshCBar = vp.plot_mesh(soil_grid, "pressure head", "", False)
actors, meshCBar2 = vp.plot_mesh_cuts(soil_grid, name, 5, False)
actors.extend([meshActor])
vp.render_window(actors, "name" , [meshCBar, meshCBar2]).Start()


# meshActor, meshCBar = vp.plot_mesh(soil_grid, "water content", "", False)
# vp.render_window([meshActor], "water content", meshCBar).Start()
# vp.render_window([rootActor], name, rootCBar).Start()
