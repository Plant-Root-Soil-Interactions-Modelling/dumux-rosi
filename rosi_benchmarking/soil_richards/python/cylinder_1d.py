# Benchmark 3 (in 1D)
#
# compares the dumux solution to the analytical solution (Figure 4abc Vanderborght et al 2005)
#
# D. Leitner, 2018
#

import os
import matplotlib.pyplot as plt
#from analytic_b3 import *
from vtk_tools import *
import van_genuchten as vg

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/rosi_benchmarking/soil_richards")

# run dumux
os.system("./richards1d input/cylinder_1d.input")


#ex = []  # list for data for export


for i in range(0, 3):
    s_, p_, z_ = read1D_vtp_data("cylinder_1d-0000" + str(i + 1) + ".vtp")
    h_ = vg.pa2head(p_)
    theta_ = vg.water_content(h_, sand)
    #ax1.plot(theta_, z_ * 100, "r+")
    #ex.append(z_)
    #ex.append(theta_)


#np.savetxt("cylinder_1d", np.vstack(ex), delimiter = ",")

plt.show()

