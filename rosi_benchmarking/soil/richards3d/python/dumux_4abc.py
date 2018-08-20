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
# os.chdir( path )
# os.system( "mpirun -n 8 ./richards3d input/b3a.input")
# os.system( "mpirun -n 8 ./richards3d_ug input/b3a_ug.input")
# 
# os.system( "mpirun -n 8 ./richards3d input/b3b.input")
# os.system( "mpirun -n 8 ./richards3d_ug input/b3b_ug.input")
# 
# os.system( "mpirun -n 8 ./richards3d input/b3c.input")
# os.system( "mpirun -n 8 ./richards3d_ug input/b3c_ug.input")

# result dumux jan3a (Figure 4a)
for i in range(0,3):
    s_, p_, z_ = read3Dp_vtp_data(path+"s0008-p000","-b3a_ug-0000"+str(i+1), 8)   
    z_ = z_*100 - 200  # m -> cm,  - offset 200 cm
    h_ = vg.pa2head(p_)
    theta_ = vg.water_content(h_, sand)
    ax1.plot(theta_,z_, "r+")   
  
# # result dumux jan3b (Figure 4b)
# for i in range(0,1):
#     s_, p_, z_ = read3Dp_vtp_data(path+"s0008-p000","-b3b_ug-0000"+str(i+1), 8) 
#     z_ = z_*100 - 200  # m -> cm,  - offset 200 cm    
#     h_ = vg.pa2head(p_) 
#     theta_ = vg.water_content(h_, loam)
#     ax2.plot(theta_,z_, "r+")
     
# result dumux jan3c (Figure 4c)
for i in range(0,1):
    s_, p_, z_ = read3Dp_vtp_data(path+"s0008-p000","-b3c_ug-0000"+str(i+1), 8) 
    z_ = z_*100 - 200  # m -> cm,  - offset 200 cm    
    h_ = vg.pa2head(p_) 
    theta_ = vg.water_content(h_, clay)
    ax3.plot(theta_,z_, "r+")

plt.show()


# os.system( "mpirun -n 8 ./richards3d_ug input/b2_ug.input")  
# s_, p_, z_  = read3Dp_vtp_data(path+"s0008-p000","-b2_ug-00001", 8)

