#
# compares the dumux solution to the analytical solution (Figure 2abc Vanderborght et al 2005)
#
# D. Leitner, 2018
#

import os
import matplotlib.pyplot as plt
# from analytic_5abcd import *
from vtk_tools import *
from van_genuchten import *
import numpy as np

sand = Parameters(0.045, 0.43, 0.15, 3, 1.1574e-04*100*3600*24)
loam = Parameters(0.08, 0.43, 0.04, 1.6, 5.7870e-06*100*3600*24)
clay = Parameters(0.1, 0.4, 0.01, 1.1, 1.1574e-06*100*3600*24)

# manually set absolute path
path = "/home/daniel/workspace/DUMUX/dumux-rosi/build-cmake/rosi_benchmarking/richards1d/"

# run dumux 
os.chdir( path )
# os.system( "./richards1d input/b4a.input")
# os.system( "./richards1d input/b4b.input")
# os.system( "./richards1d input/b4c.input")
# os.system( "./richards1d input/b4d.input")

file = open("benchmark4a.txt",'r')
d = np.loadtxt(file, delimiter=',')
t1_ = d[:,0]/(24*3600)
h1_= pa2head(d[:,1])
file.close()

file = open("benchmark4b.txt",'r')
d = np.loadtxt(file, delimiter=',')
t2_ = d[:,0]/(24*3600)
h2_= pa2head(d[:,1])
file.close()
soil = loam

file = open("benchmark4c.txt",'r')
d = np.loadtxt(file, delimiter=',')
t3_ = d[:,0]/(24*3600)
h3_= pa2head(d[:,1])
file.close()
soil = loam

file = open("benchmark4d.txt",'r')
d = np.loadtxt(file, delimiter=',')
t4_ = d[:,0]/(24*3600)
h4_= pa2head(d[:,1])
file.close()
soil = clay

soil = [sand, loam, loam, clay]

rho = 1000
g = 9.81

Eact1 = np.zeros(len(h1_))
for i,h in enumerate(h1_):    
     Eact1[i] = min(0.1,hydraulic_conductivity(h,soil[0]) * (-h + 1))

Eact2 = np.zeros(len(h2_))
for i,h in enumerate(h2_):    
    Eact2[i] = min(0.1,hydraulic_conductivity(h,soil[1]) * (-h + 1))

Eact3 = np.zeros(len(h3_))
for i,h in enumerate(h3_):    
    Eact3[i] = min(0.3,hydraulic_conductivity(h,soil[2]) * (-h + 1)) 

Eact4 = np.zeros(len(h4_))
for i,h in enumerate(h4_):    
    Eact4[i] = hydraulic_conductivity(h,soil[3]) * (-h - 1) 

# ax1.plot(t1_,Eact1,"r:")
# ax2.plot(t2_,Eact2,"r:")
# ax3.plot(t3_,Eact3,"r:")
# ax4.plot(t4_,Eact4,"r:")


plt.plot(t4_,Eact4,"r:")
plt.show()


# os.system( "./richards1d input/b4b.input")
# os.system( "./richards1d input/b4c.input")
# os.system( "./richards1d input/b4d.input")


# # Figure 5a
# s_, p_ = read1D_vtp_data(path+"benchmark4a-0000"+str(i+1)+".vtp", False)
# z_ = np.linspace(0,-200,len(s_))
# h_ = vg.pa2head(p_) 
# theta_ = vg.water_content(h_, sand)
# ax1.plot(theta_,z_, "r+")   
#     
# 
# 
# plt.show()

