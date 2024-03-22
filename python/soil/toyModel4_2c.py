import sys; sys.path.append("../modules/"); sys.path.append("../../../CPlantBox/");  sys.path.append("../../build-cmake/cpp/python_binding/")
sys.path.append("../../../CPlantBox/src/python_modules")


from rosi_richardsnc import RichardsNCSP  # C++ part (Dumux binding)
#from rosi_richards10c import Richards10CSP  # C++ part (Dumux binding)

from richards import RichardsWrapper  # Python part

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()

def write_file_array(name, data, space =","):
    name2 = './results/'+ name+ '.txt'
    with open(name2, 'a') as log:
        log.write(space.join([num for num in map(str, data)])  +'\n')
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
usemoles = True
ncsp = True
if ncsp:
    s = RichardsWrapper(RichardsNCSP())
else:
    s = RichardsWrapper(Richards10CSP())

s.initialize()

# theta_r, theta_s, alpha, n, Ks
loam = [0.045, 0.43, 0.04, 1.6, 50]

min_b = [-5., -5, -20.] 
max_b = [5., 5, 0.]
cell_number = [5, 5, 20]
# s.createGrid([0.02], [0.6], [nCells])  # [cm]
s.createGrid(min_b, max_b, cell_number, False)
cell_number = np.array(cell_number)
s.setParameter( "Soil.Grid.Cells", str(cell_number)[1:(len(str(cell_number))-1)])

s.setHomogeneousIC(-100.)  # cm pressure head
s.setOuterBC("constantFlux", 0)  #  [cm/day]
s.setInnerBC("constantFlux", -10)  #  [cm/day]
s.setParameter("Problem.doSoluteFlow", "1")

#s.setICZ_solute(0.)  # [kg/m2] 

molarMassWat = 18. # [g/mol]
densityWat = 1. #[g/cm3]
# [mol/cm3] = [g/cm3] /  [g/mol] 
molarDensityWat =  densityWat / molarMassWat # [mol/cm3] 

solidDensity = 2700 # [kg/m^3 solid]
solidMolarMass = 60.08e-3 # [kg/mol] 
# [mol / m3 solid] =[kg/m^3 solid] / [kg/mol] 
solidMolDensity = solidDensity/solidMolarMass
# [mol / m3 scv] = [mol / m3 solid] * [m3 solid /m3 space]
bulkDensity_m3 = solidMolDensity*(1.-0.43)

MolarMass = 1.8e-2 #[kg/mol] 0.02003 #[kg/mol]
exud = 1.*0# [mol/cm2/d]#1.0* MolarMass *1000# [mol/cm2/d] * [kg/mol] * [g/kg] =  [g/cm2/d]
Ds = 1e-8*0 # m^2/s
Dl = 1e-9*0

numComp = 8
numFluidComp = 2

solidDensity = 2700 # [kg/m^3 solid]
solidMolarMass = 60.08e-3 # [kg/mol] 
# [mol / m3 solid] =[kg/m^3 solid] / [kg/mol] 
solidMolDensity = solidDensity/solidMolarMass
# [mol / m3 scv] = [mol / m3 solid] * [m3 solid /m3 space]
bulkDensity_m3 = solidMolDensity*(1.-0.43)


s.setVGParameters([loam])

s.setParameter( "Soil.MolarMass", str(solidMolarMass))
s.setParameter( "Soil.solidDensity", str(solidDensity))
s.setParameter("Component.MolarMass", str(1.8e-2)) #m^2/s


s.setParameter("Problem.EnableGravity", "false")
s.setParameter("Problem.verbose", "1")
s.setParameter("Flux.UpwindWeight", "0.5")
s.setParameter("Newton.EnableAbsoluteResidualCriterion", "true")
s.setParameter("Newton.MaxAbsoluteResidual", "1.e-10")
s.setParameter("Newton.EnableChop", "true")
s.setParameter("Newton.EnableResidualCriterion", "true")
s.setParameter("Newton.EnableShiftCriterion", "true")
s.setParameter("Newton.MaxAbsoluteResidual", "1e-10")

s.setParameter("Newton.MaxRelativeShift", "1e-10")

s.setParameter("Newton.MaxSteps", "30")
s.setParameter("Newton.ResidualReduction", "1e-10")
s.setParameter("Newton.SatisfyResidualAndShiftCriterion", "true")
s.setParameter("Newton.TargetSteps", "10")
s.setParameter("Newton.UseLineSearch", "false")
s.setParameter("Newton.EnablePartialReassembly", "true")
s.setParameter("Grid.Overlap", "0")  #no effec5

s.initializeProblem()
s.setParameter("Flux.UpwindWeight", "1")
s.setCriticalPressure(-15000)  # cm pressure head


times = [0., 5./24./3600, 10./24./3600]  # days
s.ddt = 1.e-5

points = s.getCellCenters()#.flatten()

x = np.array(s.getSolutionHead())
write_file_array("pressureHead",x.flatten())
write_file_array("coord",points)

for i, dt in enumerate(np.diff(times)):

    if rank == 0:
        print("*****", "external time step", dt, " d, simulation time", s.simTime, "d, internal time step", s.ddt, "d")

    s.solve(dt)

    x = np.array(s.getSolutionHead())
    write_file_array("pressureHead",x.flatten())
    #print(x.flatten())
    write_file_array("coord",points)


