#
# compares the dumux solution to the analytical solution (Figure 2abc Vanderborght et al 2005)
#
# D. Leitner, 2018
#

import os
import matplotlib.pyplot as plt
from analytic_4abc import *
from vtk_tools import *
import van_genuchten as vg

# manually set absolute path
path = "/home/daniel/workspace/DUMUX/dumux-rosi/build-cmake/rosi_benchmarking/richards3d/"

# run dumux 
os.chdir( path )
os.system( "./richards3d input/jan3a.input")
os.system( "./richards3d input/jan3b.input")
os.system( "./richards3d input/jan3c.input")

# result dumux jan3a (Figure 4a)
for i in range(0,1):
    s_, p_, z_ = read3D_vtp_data(path+"jan3a-0000"+str(i)+".vtu", False)
    z_ = z_*100 - 200  # m -> cm,  - offset 200 cm
    h_ = vg.pa2head(p_)
    theta_ = vg.water_content(h_, sand)
    ax1.plot(theta_,z_, "r+")   
  
# result dumux jan3b (Figure 4b)
for i in range(0,3):
    s_, p_, z_ = read3D_vtp_data(path+"jan3b-0000"+str(i+1)+".vtu", False)
    z_ = z_*100 - 200  # m -> cm,  - offset 200 cm    
    h_ = vg.pa2head(p_) 
    theta_ = vg.water_content(h_, loam)
    ax2.plot(theta_,z_, "r+")
   
# # result dumux jan3c (Figure 4c)
for i in range(0,3):
    s_, p_, z_ = read3D_vtp_data(path+"jan3c-0000"+str(i+1)+".vtu", False)
    z_ = z_*100 - 200  # m -> cm,  - offset 200 cm    
    h_ = vg.pa2head(p_) 
    theta_ = vg.water_content(h_, clay)
    ax3.plot(theta_,z_, "r+")

plt.show()

