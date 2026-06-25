import sys; 
sys.path.append("../../experimental/fixedPointIter2/modules"); 
sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

from rosi_richardsnc import RichardsNCSP as RichardsNCSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part

import matplotlib.pyplot as plt
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()
import numpy as np
import functional.van_genuchten as vg

""" 
How to check the soil mass and volume balance for each cell
"""

def setSoilParam(s):    
    """ save the soil parameters
        @param: the dumux soil object
    """
    s.solidDensity = 2650 #[kg/m^3 solid] 
    s.solidMolarMass = 60.08e-3 # [kg/mol] 
    loam = [0.08, 0.43, 0.04, 1.6, 5]  # K = 5 !
    s.soil =  loam 
    s.vg_soil = vg.Parameters(s.soil) 
    # [mol / m3 solid] =[kg/m^3 solid] / [kg/mol] 
    s.solidMolDensity = s.solidDensity/s.solidMolarMass
    # [mol / m3 scv] = [mol / m3 solid] * [m3 solid /m3 space]
    s.bulkDensity_m3 = s.solidMolDensity*(1.- s.vg_soil.theta_S)
    s.bulkMassDensity_gpercm3 = s.solidDensity*(1.- s.vg_soil.theta_S)*1000/1e6

    s.setParameter( "Soil.MolarMass", str(s.solidMolarMass))
    s.setParameter( "Soil.solidDensity", str(s.solidDensity))
    s.setParameter( "Soil.css1Function", str(5))
    #decay
    s.setParameter("Soil.vmax_decay", str(0)) #mol C / m^3 scv / s 
    s.setParameter("Soil.km_decay", str(1)) #mol C / m^3 scv
   

    #sorption
    s.setParameter("Soil.CSSmax", str(0)) #[mol/cm3] 
    s.setParameter("Soil.kads", str(1)) #[cm3/mol/d]
    s.setParameter("Soil.kdes", str(1))  #[1/d]
    
    s.setVGParameters([s.soil])
    
    return s

def solve(simtimes):

    s = RichardsWrapper(RichardsNCSP())
    s = setSoilParam(s)
    s.initialize()
    s.createGrid([-5., -5., -10.], [5., 5., 0.], [2,2,2], periodic = True)  # [cm]
    #s.setVGParameters([loam])

    # theta = 0.378, benchmark is set be nearly fully saturated, so we don't care too much about the specific values
    s.setHomogeneousIC(-5000., True)  # cm pressure head

    #s.setTopBC("constantFlux", -2)  #  [cm/day] "noFlux")#
    evap = -0.1 
    s.setTopBC("atmospheric", 0.5, [[0.0, 1.0e10], [evap, evap]]) 
    s.setBotBC("constantFlux", 0.) # "noFlux")#
    # s.setParameter("Soil.BC.Top.SType", "2")  # michaelisMenten=8 (SType = Solute Type)
    # s.setParameter("Soil.BC.Top.CValue", "1.e-4")  # michaelisMenten=8 (SType = Solute Type)
    # s.setParameter("Soil.BC.Bot.SType", "6")  # michaelisMenten=8 (SType = Solute Type)
    # s.setParameter("Soil.BC.Bot.CValue", "0.")

    s.setParameter("Soil.IC.C", "0")  # g / cm3  # TODO specialised setter?

    s.setParameter("Component.MolarMass", "1.8e-2") 
    s.setParameter("Soil.MolarMass", "1.8e-2") 

    s.setParameter("Component.LiquidDiffusionCoefficient", "1.e-9")  # m2 s-1
    s.initializeProblem(maxDt = 0.01)

    # dummy sources of water and solutes to test the mass balance
    source_map = { 0: 0.01, 1:-0.01, 2: 0.02, 3:-0.02, 4: 0.03, 5:-0.03}
    #s.setSource(source_map)
    #source_map = { 0: 0.01, 1:-0.01, 2: 0.02, 3:-0.02, 4: 0.03, 5:-0.03}
    #s.setSource(source_map, eq_idx = 1)
    
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
        print('getSolutionHead',s.getSolutionHead())
        s.solve(dt, saveInnerFluxes_ = True)        
        print('getSolutionHead',s.getSolutionHead())
        Wvolafter  = s.getCellVolumes() * s.getWaterContent() # cm3
        Smassafter = s.getSolution(1)   * Wvolafter # kg
        
            
        scvIntegratedFlows   = s.getFluxesPerCell(0) * dt # cm3
        scvIntegratedFlowsS  = s.getFluxesPerCell(1) * dt # g
        scvSources  = s.getSource(0) * s.getCellVolumes() * dt # cm3
        scvSourcesS = s.getSource(1) * s.getCellVolumes() * dt # g
        
        if rank == 0:
            print('\tWvolbefore:',Wvolbefore)
            print('\tChange in water volume [cm3] per voxel:',Wvolafter-Wvolbefore)
            print('\tscvIntegratedFlows',scvIntegratedFlows)
            print('\tRMSE in water volume [cm3] per voxel:',scvIntegratedFlows+Wvolafter-Wvolbefore-scvSources  )
            # print('\tChange in solute mass [g] per voxel:',Smassafter-Smassbefore)
            print('\n\n\tRMSE for water volume balance [cm3]:',np.mean(np.abs(scvIntegratedFlows+(Wvolafter-Wvolbefore)-scvSources)), 
                            sum(scvIntegratedFlows),sum(Wvolafter-Wvolbefore),sum(scvSources))



if __name__ == "__main__":


    simTimes = [0.5]#,0.78,1]  # days
    solve(simTimes)

