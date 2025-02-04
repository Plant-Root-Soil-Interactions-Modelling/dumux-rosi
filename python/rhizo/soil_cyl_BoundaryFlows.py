import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

from rosi_richardsnc_cyl import RichardsNCCylFoam  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part

import matplotlib.pyplot as plt
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()
import numpy as np

""" 
Get the inner and outer boundary flows for the 1d axisymmetric soil
"""

r_in = 2  # cm
r_out = 10
length = 3

def solve(simtimes, N):

    loam = [0.08, 0.43, 0.04, 1.6, 5]  # K = 5 !
    s = RichardsWrapper(RichardsNCCylFoam())

    s.initialize()
    s.createGrid1d(np.linspace(r_in, r_out, N))  # [cm]
    s.setVGParameters([loam])

    # theta = 0.378, benchmark is set be nearly fully saturated, so we don't care too much about the specific values
    s.setHomogeneousIC(-1000.)  # cm pressure head

    s.setTopBC("constantFluxCyl",0.1)  #  [cm/day] "noFlux")#
    s.setBotBC("constantFluxCyl", -0.2) # "noFlux")#
    s.setParameter("Soil.BC.Top.SType", "3")  # michaelisMenten=8 (SType = Solute Type)
    s.setParameter("Soil.BC.Top.CValue", "-0.05")  # michaelisMenten=8 (SType = Solute Type)
    s.setParameter("Soil.BC.Bot.SType", "3")  # michaelisMenten=8 (SType = Solute Type)
    s.setParameter("Soil.BC.Bot.CValue", "0.1")

    s.setParameter("Soil.IC.C", "0")  # g / cm3  # TODO specialised setter?

    s.setParameter("Component.MolarMass", "1.8e-2") 

    s.setParameter("Component.LiquidDiffusionCoefficient", "1.e-9")  # m2 s-1
    s.initializeProblem(maxDt = 0.01)

    cellVolumes = s.getCellSurfacesCyl() * length # cm3
    if rank == 0:
        print(s)

    s.ddt = 1.e-3  # days

    simtimes.insert(0, 0)
    dt_ = np.diff(simtimes)
    
    for r, dt in enumerate(dt_):
        

        time = simtimes[r] 
        print('time',time)
            
        #if time >= 5:
        #    s.setSoluteTopBC([1], [0.])

        if rank == 0:
            print("*****", "#", r, "external time step", dt, " d, simulation time", s.simTime, "d, internal time step", s.ddt, "d")

        Wvolbefore = cellVolumes * s.getWaterContent() # cm3
        Smassbefore = s.getSolution(1) * (Wvolbefore / 1.e3 ) # kg
        
        s.solve(dt)
        
        Wvolafter = cellVolumes*s.getWaterContent() # cm3
        Smassafter = s.getSolution(1) * (Wvolafter / 1.e3 ) # kg
        
            
        rootSoilFluxes = s.getInnerFlow(0, length) * dt # cm3
        rootSoilFluxesS = s.getInnerFlow(1, length) * dt # kg
        soilSoilFluxes = s.getOuterFlow(0, length) * dt # cm3
        soilSoilFluxesS = s.getOuterFlow(1, length) * dt # kg
        
        # TODO: currently, setSource not properly implemented for richards and richards 2c
        # so left out of the mass balance.
        # scvSources = s.getSource(0) * cellVolumes * dt # cm3
        # scvSourcesS = s.getSource(1) * cellVolumes * dt # kg
        
        if rank == 0:
            print('\tChange in water volume [cm3] per voxel:',sum(Wvolafter-Wvolbefore))
            print('\tChange in solute mass [kg] per voxel:',sum(Smassafter-Smassbefore))
            print('\n\n\tRMSE for water volume balance [cm3]:',np.mean(np.abs(rootSoilFluxes+soilSoilFluxes+sum(Wvolafter-Wvolbefore))))
            print('\tRMSE for solute mass balance [kg]:',np.mean(np.abs(rootSoilFluxesS+soilSoilFluxesS+sum(Smassafter-Smassbefore))),'\n\n')
            

       


if __name__ == "__main__":


    simTimes = [0.5,0.78,1]  # days
    solve(simTimes, 2)

