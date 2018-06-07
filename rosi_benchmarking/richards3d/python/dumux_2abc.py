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

# run dumux
os.chdir( path )
# os.system( "./richards3d_steadystate input/b1a.input")
# os.system( "./richards3d_steadystate input/b1b.input")
# os.system( "./richards3d_steadystate input/b1c.input")
# os.system( "./richards3d_steadystate input/b1a_ug.input")
# os.system( "./richards3d_steadystate input/b1b_ug.input")
# os.system( "./richards3d_steadystate input/b1c_ug.input")
# os.system( "./richards3d_steadystate input/b1a_ug2.input")
# os.system( "./richards3d_steadystate input/b1b_ug2.input")
# os.system( "./richards3d_steadystate input/b1c_ug2.input")

j = 2

if j==0:
    print("high resolution, regular")
if j==1:
    print("unstructured (distmesh)")
if j==2:
    print("unstructured with hr layer (distmesh)")


# Result dumux (Figure 2a)
if j==0:
    s_, p_, z_ = read3D_vtp_data(path+"b1a-00001.vtu", False)
if j==1:
    s_, p_, z_ = read3D_vtp_data(path+"b1a_ug-00001.vtu", False)
if j==2:

    s_, p_, z_ = read3D_vtp_data(path+"b1a_ug2-00001.vtu", False)

z_ = z_*100 - 200  # m -> cm,  - offset 200 cm
h_ = vg.pa2head(p_) 
ax1.plot(h_,z_, "r+")
print("dof : ", len(z_)) # in the parallel code this value is not correct
 
# result dumux jan1 (Figure 2b)
if j==0:
    s_, p_, z_ = read3D_vtp_data(path+"b1b-00001.vtu", False)
if j==1:
    p_, z_ = np.array([0]), np.array([0])
    pass
    # s_, p_, z_ = read3D_vtp_data(path+"b1b_ug-00001.vtu", False)
if j==2:
    p_, z_ = np.array([0]), np.array([0])  
    pass    
    # s_, p_, z_ = read3D_vtp_data(path+"b1b_ug2-00001.vtu", False)

z_ = z_*100 - 200  # m -> cm,  - offset 200 cm
h_ = vg.pa2head(p_) 
ax2.plot(h_,z_, "r+")
print("dof : ", len(z_)) # in the parallel code this value is not correct
 
# result dumux jan1 (Figure 2c)
if j==0:
    s_, p_, z_ = read3D_vtp_data(path+"b1c-00001.vtu", False)
if j==1:
    s_, p_, z_ = read3D_vtp_data(path+"b1c_ug-00001.vtu", False)
if j==2:
    s_, p_, z_ = read3D_vtp_data(path+"b1c_ug2-00001.vtu", False)

z_ = z_*100 - 200  # m -> cm,  - offset 200 cm
h_ = vg.pa2head(p_) 
ax3.plot(h_,z_, "r+")
print("dof : ", len(z_)) # in the parallel code this value is not correct

plt.show()

