import sys; sys.path.append("../modules/"); sys.path.append("../../../CPlantBox/");  sys.path.append("../../build-cmake/cpp/python_binding/")
sys.path.append("../../../CPlantBox/src/python_modules")

from rosi_richards_cyl import RichardsCylFoam  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part

import matplotlib.pyplot as plt
import numpy as np
import os

def write_file_array(name, data):
    name2 = './results/'+ name+ '.txt'
    with open(name2, 'a') as log:
        log.write(','.join([num for num in map(str, data)])  +'\n')

""" 
Cylindrical 1D model, diffusion only (DuMux), Michaelis Menten

everything scripted, no input file needed, also works parallel with mpiexec
"""

results_dir="./results/"
if not os.path.exists(results_dir):
    os.makedirs(results_dir)
else:
    test = os.listdir(results_dir)
    for item in test:
        try:
            os.remove(results_dir+item)
        except:
            pass


s = RichardsWrapper(RichardsCylFoam())
s.initialize()

#theta_r, theta_s, alpha, n, Ks?
loam = [0.045, 0.43, 0.04, 1.6, 50]

#s.createGrid([0.02], [0.6], [500])  # [cm]
# with open("Hp_results.txt") as f:
    # list2 = [row.split(',')[0] for row in f]
# list2 = list2[1:]
# list2 = [float(ll) for ll in list2]

# s.createGrid1d(np.array(list2) + 0.02)


s.createGrid([0.02], [0.6], [700]) 
s.setHomogeneousIC(-100.)  # cm pressure head
s.setOuterBC("fluxCyl", 0.)  #  [cm/day]
s.setInnerBC("constantFluxCyl", -1)  #  [cm/day]

s.setVGParameters([loam])
s.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
s.setParameter("Newton.MaxAbsoluteResidual", "1.e-10")
s.initializeProblem()

s.setCriticalPressure(-15000)  # cm pressure head

times = [10,10]  # days
s.ddt = 1.e-5

points = s.getDofCoordinates()

x = np.array(s.getSolutionHead())
write_file_array("pressureHead",x.flatten())
write_file_array("coord",points.flatten())
write_file_array("theta",np.array(s.getWaterContent()).flatten())
write_file_array("getSaturation",np.array(s.getSaturation()).flatten())
write_file_array("krs",np.array(s.getKrw()).flatten())
#raise Exception
for i, dt in enumerate(times):

    print("*****", "external time step", dt, " d, simulation time", s.simTime, "d, internal time step", s.ddt, "d")

    s.solve(dt)

    points = s.getDofCoordinates()

    x = np.array(s.getSolutionHead())
    write_file_array("pressureHead",x.flatten())
    write_file_array("coord",points.flatten())
    write_file_array("theta",np.array(s.getWaterContent()).flatten())
    write_file_array("getSaturation",np.array(s.getSaturation()).flatten())
    write_file_array("krs",np.array(s.getKrw()).flatten())
    #print(points.flatten())
    #print(x.flatten())
    #print(np.array(s.getWaterContent()).flatten())


