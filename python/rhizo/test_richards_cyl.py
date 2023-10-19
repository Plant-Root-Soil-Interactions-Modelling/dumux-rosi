import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

import functional.van_genuchten as vg
from fv.fv_grid import *
import fv.fv_richards as richards  # Python solver

from rosi_richards_cyl import RichardsCylFoam  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()  # MPI

import matplotlib.pyplot as plt
import numpy as np
import os

SMALL_SIZE = 22
MEDIUM_SIZE = 22
BIGGER_SIZE = 22
plt.rc('font', size = SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize = SMALL_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize = MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize = SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize = BIGGER_SIZE)  # fontsize of the figure title

""" 
Compares Dumux and Python cylindircal richards solutions 
for our rhizosphere setting (as used in coupled_c12_rhizo.py, mode = "dumux" and "python")

Dumux performs better! Python solver fails for higher resolutions or larger time steps (sand does not converge). 
Dumux solution always converges, but can be wrong for large times steps. 
"""

sand = [0.045, 0.43, 0.15, 3, 1000]
loam = [0.08, 0.43, 0.04, 1.6, 50]
clay = [0.1, 0.4, 0.01, 1.1, 10]
soil = loam

""" fixed flux """
j = -0.1  # [cm3 / day]

# Parameters
times = np.linspace(0, 2.5, 11)  # days
x0 = -100  # cm
wilting_point = -15000

# grid
dof = 20  # better with 20
NC = dof + 1
a_in = 0.02
a_out = 0.1
lb = 1.1  # logbase
points = np.logspace(np.log(a_in) / np.log(lb), np.log(a_out) / np.log(lb), NC, base = lb)
# points = np.linspace(a_in, a_out, NC)

# DUMUX solver (intialization as in rhizo.models initialize_dumux_)
s = RichardsWrapper(RichardsCylFoam())
s.initialize()
s.createGrid1d(points)
s.setHomogeneousIC(x0)  # cm pressure head
s.setVGParameters([soil])
s.setOuterBC("fluxCyl", 0.)  #  [cm/day]
s.setInnerBC("fluxCyl", j / (2 * np.pi * a_in * 1.))  # [cm3 / day] -> [cm3 / cm3 /day]
s.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
s.setParameter("Newton.MaxAbsoluteResidual", "1.e-10")
s.setParameter("TimeLoop.MaxTimeStepSize", "60")  # seconds
s.initializeProblem()
s.setCriticalPressure(wilting_point)  # cm pressure head

# PYTHON solver (intialization as in rhizo.models initialize_python_)
s2 = richards.FVRichards1D(FVGrid1Dcyl(points), soil)
s2.x0 = np.ones((NC - 1,)) * x0
s2.bc[(0, 0)] = ("flux_in_out", [j / (2 * np.pi * a_in * 1.), wilting_point])  # # [cm3 / day] -> [cm3 / cm3 /day]

dt_ = np.diff(times)
col = ["--r*", "--b*", "--g*", "--c*", "--m*", "--y*", ] * 2
col2 = [":r*", ":b*", ":g*", ":c*", ":m*", ":y*", ] * 2
for i, dt in enumerate(dt_):
    points = s.getDofCoordinates()
    x = s.getSolutionHead()
    points2 = s2.grid.centers()
    x2 = s2.x0
    plt.plot(points[:], x, col[i % len(col)], label = "dumux {} days".format(s.simTime))
    plt.plot(points2[:], x2, col2[i % len(col2)], label = "python {} days".format(s.simTime))
    s.solve(dt)
    s2.solve([dt], 60 / 3600 / 24)  # [s] -> [day]

plt.xlabel("distance from root axis (cm)")
plt.xlabel("soil matric potential (cm)")
# plt.legend()
plt.show()

