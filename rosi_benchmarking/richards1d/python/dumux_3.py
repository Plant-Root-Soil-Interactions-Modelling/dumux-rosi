#
# compares the dumux solution to the analytical solution (Figure 3, Vanderborght et al 2005)
#
# D. Leitner, 2018
#
import os
import matplotlib.pyplot as plt
from vtk_tools import *
import van_genuchten as vg
from scipy import integrate

loam = vg.Parameters(0.08, 0.43, 0.04, 1.6, 5.7870e-06*100*3600*24) # Ksat = 50

Jw = -0.5 # cm/day

K = lambda psi: vg.hydraulic_conductivity(psi,loam)
F = lambda psi: 1./(Jw/K(psi) - 1.)

psi_ = np.linspace(0,-300,300) # psi(-54) = 0
dz = np.zeros(len(psi_),)
for i in range(0, len(psi_)):
    ans,err = integrate.quad(F,0,psi_[i])
    dz[i] = ans

z1 = dz + (-53); # this value is not clear to me any more 
plt.plot(psi_, z1)

plt.xlabel('\psi (cm)');
plt.ylabel('Depth (cm)');
plt.xlim(-300,0)
plt.ylim(-60,0)

# manually set absolute path
path = "/home/daniel/workspace/DUMUX/dumux-rosi/build-cmake/rosi_benchmarking/richards1d/"

# run dumux 
os.chdir(path)
os.system( "./richards1d input/jan2.input")

# result dumux jan1 (Figure 2a)
s_, p_ = read1D_vtp_data(path+"jan2-00001.vtp", False)
z_ = np.linspace(0,-54,len(s_))
h_ = vg.pa2head(p_) 
plt.plot(h_,z_, "r+")
 
plt.show()
