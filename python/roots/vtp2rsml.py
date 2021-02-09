import sys; sys.path.append("../modules/"); sys.path.append("../../../CPlantBox/");  sys.path.append("../../../CPlantBox/src/python_modules")

from xylem_flux import XylemFluxPython  # Python hybrid solver
import plantbox as pb
import rsml_reader as rsml
import rsml_writer as rsmlw
import vtk_tools as vt
import vtk_plot as vp

from math import *
import numpy as np
import matplotlib.pyplot as plt

""" 
Converts a DuMux output vtp to a RSML
"""

file_in = "../../grids/RootSystem8.vtp"
file_out = "../../grids/RootSystem8.rsml"

""" read vtp """
pd = vt.read_vtp(file_in)

""" node data m -> cm """
nodes = vt.np_points(pd)
pd.SetPoints(vt.vtk_points(100.*nodes))

""" convert from cell to point data """
# print(pd.GetCellData())
n = pd.GetPointData().GetNumberOfArrays()
for i in range(0,n):
    pd.GetPointData().RemoveArray(0)
    
segs = vt.np_cells(pd)

# age (age [s] -> emergence_time [day])
age_id = 6
age = np.zeros((segs.shape[0] + 1,))
age[1:], _ = vt.np_data(pd, age_id, True)
age[0] = age[1]
max_age = np.max(age)
et = max_age - age
pd.GetPointData().AddArray(vt.vtk_data(et))

# radius (radius [m] -> diameter [cm])
radius_id = 3
radii = np.zeros((segs.shape[0] + 1,))
radii[1:], _ = vt.np_data(pd, radius_id, True)
radii[0] = radii[1]
pd.GetPointData().AddArray(vt.vtk_data(2. * radii * 100))

# type (order = type)
order_id = 4
types = np.zeros((segs.shape[0] + 1,))
types[1:], _ = vt.np_data(pd, order_id, True)
types[0] = types[1]
types = (types > 0) + np.ones(types.shape)  # <---------------------------
pd.GetPointData().AddArray(vt.vtk_data(types))

meta = rsmlw.Metadata()
meta.set_fun_names(["emergence_time", "diameter", "type"])

vt.write_rsml(file_out, pd, meta)

print("fin")
