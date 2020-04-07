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
#os.chdir("../../../build-cmake/rosi_benchmarking/soil_tracer")

# run dumux
#os.system("./tracer_cyl input/cylinder_1d.input")

#read in data
os.chdir("../../../rosi_benchmarking/soil_tracer/python")
data = np.loadtxt("cylinder_1d_Comsol.txt", skiprows=8)
z_comsol = data[:,0]
c_comsol = data[:,-1]

# Figure 
#s_, c_, z_ = read1D_vtp_data("cylinder_1d-00001.vtp", 13)
#plt.plot(c_, z_ * 100, "r:")
plt.plot(z_comsol, c_comsol, "r:")

plt.xlabel('distance from the root surface (cm)')
plt.ylabel('phosphate concentration (mol/cm³)')
#plt.legend(["dumux", "comsol"], loc='lower right')
plt.show()
