import sys; sys.path.append("../../../python/modules/"); sys.path.append("../../../../CPlantBox/")

import vtk_plot as vp
import vtk_tools as vt
import rsml_writer

"""
plot DuMux .vtp output, converts it to a rsml file

TODO radius is not working, needed as PointData, called "radius"
TODO properties and functions need to be written (remove META argument, improve example) 
TODO testing 
"""

name = "soybean_Honly-00001"
# name = "dumux_c12_2cm"  # to not convert m->cm for this file

pd = vp.read_vtp(name + ".vtp")  
print(pd.GetBounds()) # xmin, xmax, ymin, ymax, zmin, zmax

# Convert from m to cm
np_points = vt.np_points(pd)
points = vt.vtk_points(np_points*100) # m -> cm
pd.SetPoints(points)

vp.plot_roots(pd, "p xylem [cm]")
#vp.plot_roots(pd, "subType")


# write RSML - this will only work for non-periodic vtp
# roots of periodic vtp will end at domain boundary (since no connection is defined in the vtp)
meta = rsml_writer.Metadata()
vt.write_rsml(name+".rsml",pd, meta, 6) # 6 is the data index of root order 

print("fin")