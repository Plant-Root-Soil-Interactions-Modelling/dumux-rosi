#
# compares the dumux solution 3d to the analytical solution (Figure 2abc Vanderborght et al 2005)
#
# D. Leitner, 2018
#

import os
import matplotlib.pyplot as plt
from analytic_2abc import *
from vtk_tools import *
import van_genuchten as vg

# manually set absolute path
path = "/home/daniel/workspace/DUMUX/dumux-rosi/build-cmake/rosi_benchmarking/richards3d/"

# # # run dumux 
os.chdir( path )
os.system( "./richards3d input/jan1a.input")
# os.system( "./richards3d input/jan1b.input")
# os.system( "./richards3d input/jan1c.input")

# result dumux jan1 (Figure 2a)
s_, p_, z_ = read3D_vtp_data(path+"jan1a-00000.vtu", False)
z_ = z_*100 - 200  # m -> cm,  - offset 200 cm
h_ = vg.pa2head(p_) 
ax1.plot(h_,z_, "r+")
 
# # result dumux jan1 (Figure 2b)
# s_, p_, z_  = read3D_vtp_data(path+"jan1b-00001.vtu", False)
# z_ = z_*100 - 200  # m -> cm,  - offset 200 cm
# h_ = vg.pa2head(p_) 
# ax2.plot(h_,z_, "r+")
#   
# # result dumux jan1 (Figure 2c)
# s_, p_, z_  = read3D_vtp_data(path+"jan1c-00001.vtu", False)
# z_ = z_*100 - 200  # m -> cm,  - offset 200 cm
# h_ = vg.pa2head(p_) 
# ax3.plot(h_,z_, "r+")

plt.show()

