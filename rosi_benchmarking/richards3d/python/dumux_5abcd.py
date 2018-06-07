#
# compares the dumux solution to the analytical solution (Figure 2abc Vanderborght et al 2005)
#
# D. Leitner, 2018
#

import os
import matplotlib.pyplot as plt
from analytic_5abcd import *
from vtk_tools import *
from van_genuchten import *
import numpy as np
from math import *

sand = Parameters(0.045, 0.43, 0.15, 3, 1.1574e-04*100*3600*24)
loam = Parameters(0.08, 0.43, 0.04, 1.6, 5.7870e-06*100*3600*24)
clay = Parameters(0.1, 0.4, 0.01, 1.1, 1.1574e-06*100*3600*24)

# manually set absolute path
path = "/home/daniel/workspace/DUMUX/dumux-rosi/build-cmake/rosi_benchmarking/richards1d/"

# run dumux 
os.chdir( path )
os.system( "./richards1d input/b4a_hr.input")
os.system( "./richards1d input/b4b_hr.input")
os.system( "./richards1d input/b4c_hr.input")
os.system( "./richards1d input/b4d_hr.input")
os.system( "./richards1d input/b4a.input")
os.system( "./richards1d input/b4b.input")
os.system( "./richards1d input/b4c.input")
os.system( "./richards1d input/b4d.input")

# open results
num = ['a','b','c','d', 'ahr','bhr','chr','dhr']
t = []
y = []


for n in num: 
    with open("benchmark4"+n+".txt",'r') as f: 
        d = np.loadtxt(f, delimiter=',')
    t_ = d[:,0]/(24*3600)
    f_ = d[:,1]*(24*3600)/1000*100; # from kg/(s m²) to cm/day
    t.append(t_)
    y.append(f_)

# prepare plot
axis = [ax1, ax2, ax3, ax4, ax1, ax2, ax3, ax4]
lt = ["r:","r:","r:","r:","r--","r--","r--","r--"]
for i,a in enumerate(axis):
        a.plot(t[i],y[i],lt[i])

plt.show()

