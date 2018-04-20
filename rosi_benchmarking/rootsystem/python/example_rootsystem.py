#
# compares the dumux solution of 1d root model to its analytical solution 
#
# D. Leitner, 2018
#

import os
import matplotlib.pyplot as plt
from vtk_tools import *
from math import *
import van_genuchten as vg

g = 9.81                  # gravitational acceleration (m/s^2)   
rho = 1.e3                # density of water, (kg/m^3)      
ref = 1.e5

def toPa(ph): # cm pressure head to Pascal (kg/ (m s^2))
    return ref + ph/100. * rho * g

def toHead(pa): # Pascal (kg/ (m s^2)) to cm pressure head
    return (pa-ref) * 100 / rho / g

#
# Dumux Solution: 
#     Geometry is in the DGF, radii, kr, & kz
#     Rest of parameters are in the dumux input file
#

# manually set absolute path
path = "/home/daniel/workspace/DUMUX/dumux-rosi/build-cmake/rosi_benchmarking/rootsystem/"

# run dumux 
# os.chdir( path )
# os.system( "./rootsystem input/rootsystem.input")


x, nodes = read3D_vtp_data(path+"rootsystem-00001.vtp", False)
x = toHead(x)
# plt.plot(h_,z_, "r+")
# plt.ylabel("Depth (m)")
# plt.xlabel("Xylem pressure (cm)")
# plt.show()


xmin = min(x)
xmax = max(x)
print("from ", xmin, "to", xmax, " cm pressure head")

for i,n in enumerate(nodes): 
    c = (x[i]-xmin)/(xmax-xmin)
    plt.plot([n[0], n[0]], [n[2], n[2]],"*",color = plt.cm.jet(c))    
plt.axis('equal')
plt.xlabel("cm")
plt.ylabel("cm")
plt.show()
