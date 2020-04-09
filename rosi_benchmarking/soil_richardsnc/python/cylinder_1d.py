#radially symmetric cylinder (in 1D)
#
# compares the dumux solution to the comsol solution
#
# Just a quick check that water movement including a tracer is plausible.
#
# D. Leitner, 2020
#

import os
import matplotlib.pyplot as plt
from vtk_tools import *
import van_genuchten as vg


# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/rosi_benchmarking/soil_richardsnc")

# run dumux
os.system("./richardsnc1d_cyl input/cylinder_1d.input")

#read in data
os.chdir("../../../build-cmake/rosi_benchmarking/soil_richardsnc/python")
data = np.loadtxt("cylinder_1d_Comsol.txt", skiprows=8)
z_comsol = data[:,0]
c_comsol1 = data[:,25] #10 days 
c_comsol2 = data[:,-1] #20 days


# Figure 
s_, c_, z_ = read1D_vtp_data("cylinder_1d-00001.vtp", 13)
plt.plot(c_, (z_ -0.0002)*100, "r:", z_comsol, c_comsol1, "b:")
plt.plot( z_comsol, c_comsol1, "b:")
s_, c_, z_ = read1D_vtp_data("cylinder_1d-00002.vtp", 13)
plt.plot(c_, (z_-0.0002) * 100, "g:", z_comsol, c_comsol2, "k:")


plt.xlabel('distance from the root surface (cm)')
plt.ylabel('phosphate concentration (g/cm³)')
plt.legend(["dumux, 10 days", "comsol, 10 days", 'dumux, 20 days", "comsol, 20 days"], loc='lower right')
plt.show()
