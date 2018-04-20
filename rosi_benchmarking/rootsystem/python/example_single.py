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
# Analytical solution 
# 

# Parameters
L = 0.5                # length of single straight root (m)
a = 2.e-3              # radius (m)
kz = 5.e-13            # axial conductivity (m^5 s / kg) (mal rho ergibt die alten einheiten)
kr = 2.e-9              # radial conductivity per root type (m^2 s / kg) 
p0 = toPa(-1000)        # dircichlet bc at top (ćm)
pL = toPa(-500)         # dircichlet bc at bot (ćm)
p_s = toPa(-200)        # static soil pressure (cm) 


# the solution of 
# d^2/dz^2 p_r = - c p_r + c p_s, 
# is 
# p_r = p_s + d_1 exp(sqrt(c) z ) + d_2 exp(-sqrt(c) z)
# with 
# c = 2 a pi kr/kz  
#

c = 2*a*pi*kr/kz

# BC 
# top: p_r(0) = p0
# bot: p_r(-L) = pL 
# bot: qz(L) = 0, -> d/dz p_r (L) = rho*g 

AA = np.array([[1,1], [sqrt(c)*exp(-sqrt(c)*L), -sqrt(c)*exp(sqrt(c)*L)] ]) # dirichlet top, neumann bot
bb = np.array([p0-p_s, -rho*g]) #

# AA = np.array([[1,1], [exp(-sqrt(c)*L), exp(sqrt(c)*L)] ]) # dirichlet top & bot
# bb = np.array([p0-p_s, pL-p_s]) #
  
d = np.linalg.solve(AA, bb) # compute constants d_1 and d_2 from bc

p_r = lambda z: toHead( p_s + d[0]*exp(sqrt(c)*z) + d[1]*exp(-sqrt(c)*z) )

za_ = np.linspace(0,-L,100)
pr = list(map(p_r, za_))

#
# Dumux Solution 
#     Geometry is in the DGF
#     Parameters are int the dumux input file
#

# manually set absolute path
path = "/home/daniel/workspace/DUMUX/dumux-rosi/build-cmake/rosi_benchmarking/rootsystem/"

# run dumux 
os.chdir( path )
# os.system( "./rootsystem input/singleroot_S1.input")
# p_ = read1D_vtp_data(path+"singleroot_S1-00001.vtp", False)
os.system( "./rootsystem input/singleroot_S2.input")
p_ = read1D_vtp_data(path+"singleroot_S2-00001.vtp", False)
z_ = np.linspace(0,-0.5,len(p_))
h_ = vg.pa2head(p_) 
plt.plot(h_,z_, "r+")
plt.plot(pr,za_)
plt.ylabel("Depth (m)")
plt.xlabel("Xylem pressure (cm)")
plt.show()




