import sys; sys.path.append("../../../python/modules/"); sys.path.append("../../../../CPlantBox/")

import vtk_plot as vp
import vtk_tools as vt
import rsml_writer

"""
plot DuMux .vtp output, converts it to a rsml file

TODO testing 

TODO properties and functions need to be written (remove META argument, improve example) 
TODO non-periodic .vtp do not work (propably if multiple basal roots are present)
"""

# name = "soybean_Honly-00001"
# name = "dumux_c12_2cm"  # to not convert m->cm for this file
name = "Glycine_max_154days"

pd = vp.read_vtp(name + ".vtp")
print(pd.GetBounds())  # xmin, xmax, ymin, ymax, zmin, zmax

# Convert from m to cm
np_points = vt.np_points(pd)
# np_points = np_points * 100 # m -> cm
points = vt.vtk_points(np_points)
pd.SetPoints(points)
print("Number of points opened: ", np_points.shape[0])

vp.plot_roots(pd, "subType")  # "p xylem [cm]"

# write RSML - this will only work for non-periodic vtp
# roots of periodic vtp will end at domain boundary (since no connection is defined in the vtp)
meta = rsml_writer.Metadata()
vt.write_rsml(name + ".rsml", pd, meta, 2)  # 6 is the data index of root order; 2 for Glycine_max_154days

print("fin")
