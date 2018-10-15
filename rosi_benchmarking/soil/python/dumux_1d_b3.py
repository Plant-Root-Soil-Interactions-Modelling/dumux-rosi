# Benchmark 3 (in 1D)
#
# compares the dumux solution to the analytical solution (Figure 4abc Vanderborght et al 2005)
#
# D. Leitner, 2018
#

import os
import matplotlib.pyplot as plt
from analytic_b3 import *
from vtk_tools import *
import van_genuchten as vg

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)
os.chdir("../../../build-cmake/rosi_benchmarking/soil")

# run dumux
# os.system("./richards1d benchmarks_1d/b3a.input")
# os.system("./richards1d benchmarks_1d/b3b.input")
# os.system("./richards1d benchmarks_1d/b3c.input")

ex = []

# result (Figure 4a)
for i in range(0, 3):
    s_, p_ = read1D_vtp_data("benchmark1d_3a-0000" + str(i + 1) + ".vtp", False)
    z_ = np.linspace(0, -200, len(s_))
    h_ = vg.pa2head(p_)
    theta_ = vg.water_content(h_, sand)
    ax1.plot(theta_, z_, "r+")
    ex.append(z_ / 100)
    ex.append(theta_)

# result (Figure 4b)
for i in range(0, 3):
    s_, p_ = read1D_vtp_data("benchmark1d_3b-0000" + str(i + 1) + ".vtp", False)
    z_ = np.linspace(0, -200, len(s_))
    h_ = vg.pa2head(p_)
    theta_ = vg.water_content(h_, loam)
    ax2.plot(theta_, z_, "r+")
    ex.append(z_ / 100)
    ex.append(theta_)

# result (Figure 4c)
for i in range(0, 3):
    s_, p_ = read1D_vtp_data("benchmark1d_3c-0000" + str(i + 1) + ".vtp", False)
    z_ = np.linspace(0, -200, len(s_))
    h_ = vg.pa2head(p_)
    theta_ = vg.water_content(h_, clay)
    ax3.plot(theta_, z_, "r+")
    ex.append(z_ / 100)
    ex.append(theta_)

np.savetxt("dumux1d", np.vstack(ex), delimiter = ",")

# print("Fully saturated at ", vg.head2pa([0.]), " Pa")
# print("Fully saturated at ", vg.pa2head([1.e5]), " cm")
# print("-200 cm at ", vg.head2pa([-200]), "Pa ")
# print("-10000 cm", vg.head2pa([-10000]))
#
# print(vg.pa2head([0.]), "cm")
# print(vg.water_content(0, loam))

plt.show()

