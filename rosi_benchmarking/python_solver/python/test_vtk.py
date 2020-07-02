import numpy as np
import vtk

from solver.xylem_flux import XylemFluxPython  # Python hybrid solver
import solver.plantbox as pb

import vtk_plot as vp
import vtk_tools as vt

""" Soil Grid """
min_b = [-4., -4., -15.]
max_b = [4., 4., 0.]
cell_number = [8, 8, 15]  # [32, 32, 60]  # [8, 8, 15]  # 32, 32, 60
soil_grid = vp.uniform_grid(np.array(min_b), np.array(max_b), np.array(cell_number))

soil_water_content = vt.vtk_data(np.array(range(0, np.prod(cell_number))) * 0.01)  # ordering... z, y, x
soil_water_content.SetName("water content")
soil_grid.GetCellData().AddArray(soil_water_content)
soil_pressure = vt.vtk_data(np.ones((np.prod(cell_number),)))
soil_pressure.SetName("pressure head")  # in macroscopic soil
soil_grid.GetCellData().AddArray(soil_pressure)

old_rs = XylemFluxPython("../grids/RootSystem_big.rsml")
ana = pb.SegmentAnalyser(old_rs.rs)
ana.filter("creationTime", 0., 8.)
ana.crop(pb.SDF_PlantBox(7.76, 7.76, 14.76))  # that's akward.. (but I wait for the final rsml).
ana.pack()
pd = vp.segs_to_polydata(ana)

# meshActor, meshCBar = vp.plot_mesh(soil_grid, "water content", win_title = "mesh plot", render = True)
# meshActor2, meshCBar2 = vp.plot_mesh_cuts(soil_grid, "water content", win_title = "cut mesh plot", render = True)
# rootActor, rootCBar = vp.plot_roots(pd, "creationTime", win_title = "root system plot", render = True)

meshActor2, meshCBar2 = vp.plot_mesh_cuts(soil_grid, "water content", win_title = "cut mesh plot", render = False)
rootActor, rootCBar = vp.plot_roots(pd, "creationTime", win_title = "root system plot", render = False)
meshActor2.extend([rootActor])
lut = meshActor2[0].GetMapper().GetLookupTable()  # same same
rootActor.GetMapper().SetLookupTable(lut)
vp.render_window(meshActor2, "mixed plot", meshCBar2).Start()
