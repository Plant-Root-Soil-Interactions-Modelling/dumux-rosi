import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

from rosi_richardsnc import RichardsNCSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part

import matplotlib.pyplot as plt
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()
import numpy as np

""" 
How to check the soil mass and volume balance for each cell
"""


def solve(simtimes):

    loam = [0.08, 0.43, 0.04, 1.6, 5]  # K = 5 !
    s = RichardsWrapper(RichardsNCSP())

    s.initialize()
    s.createGrid([-5., -5., -10.], [5., 5., 0.], [2,2,2], periodic = True)  # [cm]
    s.setVGParameters([loam])

    # theta = 0.378, benchmark is set be nearly fully saturated, so we don't care too much about the specific values
    s.setHomogeneousIC(-5000.)  # cm pressure head

    s.setTopBC("constantFlux", 0.2)  #  [cm/day] "noFlux")#
    s.setBotBC("freeDrainage") # "noFlux")#
    s.setParameter("Soil.BC.Top.SType", "2")  # michaelisMenten=8 (SType = Solute Type)
    s.setParameter("Soil.BC.Top.CValue", "1.e-4")  # michaelisMenten=8 (SType = Solute Type)
    s.setParameter("Soil.BC.Bot.SType", "6")  # michaelisMenten=8 (SType = Solute Type)
    s.setParameter("Soil.BC.Bot.CValue", "0.")

    s.setParameter("Soil.IC.C", "0")  # g / cm3  # TODO specialised setter?

    s.setParameter("Component.MolarMass", "1.8e-2") 

    s.setParameter("Component.LiquidDiffusionCoefficient", "1.e-9")  # m2 s-1
    s.initializeProblem(maxDt = 0.01)

    # dummy sources of water and solutes to test the mass balance
    source_map = { 0: 0.01, 1:-0.01, 2: 0.02, 3:-0.02, 4: 0.03, 5:-0.03}
    s.setSource(source_map)
    source_map = { 0: 0.01, 1:-0.01, 2: 0.02, 3:-0.02, 4: 0.03, 5:-0.03}
    s.setSource(source_map, eq_idx = 1)
    
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

        Wvolbefore  = s.getCellVolumes() * s.getWaterContent() # cm3
        Smassbefore = s.getSolution(1)   * Wvolbefore # g
        
        s.solve(dt, saveInnerFluxes_ = True)
        
        Wvolafter  = s.getCellVolumes() * s.getWaterContent() # cm3
        Smassafter = s.getSolution(1)   * Wvolafter # kg
        
            
        scvFluxes   = s.getFluxesPerCell(0) * dt # cm3
        scvFluxesS  = s.getFluxesPerCell(1) * dt # kg
        scvSources  = s.getSource(0) * s.getCellVolumes() * dt # cm3
        scvSourcesS = s.getSource(1) * s.getCellVolumes() * dt # kg
        
        if rank == 0:
            print('\tChange in water volume [cm3] per voxel:',Wvolafter-Wvolbefore)
            print('\tChange in solute mass [g] per voxel:',Smassafter-Smassbefore)
            print('\n\n\tRMSE for water volume balance [cm3]:',np.mean(np.abs(scvFluxes+(Wvolafter-Wvolbefore)-scvSources)))
            print('\tRMSE for solute mass balance [g]:',np.mean(np.abs(scvFluxesS+(Smassafter-Smassbefore)-scvSourcesS)),'\n\n')



if __name__ == "__main__":


    simTimes = [0.5,0.78,1]  # days
    solve(simTimes)

