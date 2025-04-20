import sys; sys.path.append("../modules"); sys.path.append("../../../build-cmake/cpp/python_binding/");
sys.path.append("../../../../CPlantBox");  sys.path.append("../../../../CPlantBox/src")

#from rosi_richards10c import Richards10CSPILU as RichardsNCSP 
from rosi_richards4c import Richards4CSPILU as RichardsNCSP # C++ part (Dumux binding)
#from rosi_richardsnc import RichardsNCSP as RichardsNCSP 
from richards import RichardsWrapper  # Python part

import matplotlib.pyplot as plt
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()
import numpy as np

""" 
How to check the soil mass and volume balance for each cell
"""


def solve(simtimes):

    loam = [0.041, 0.494, 0.0256, 1.49, 245]  # K = 5 !
    s = RichardsWrapper(RichardsNCSP(), True)

    s.initialize()
    s.createGrid([-2., -1., -1.], [0., 0., 0.], [2,1,1], periodic = True)  # [cm]
    s.setVGParameters([loam])

    # theta = 0.378, benchmark is set be nearly fully saturated, so we don't care too much about the specific values
    s.setHomogeneousIC(-1000.)  # cm pressure head
    s.setParameter( "Soil.MolarMass", str(0.06008))#0.06008 # 2650
    s.setParameter( "Soil.solidDensity", str(2650))
    s.setTopBC("constantFlux", 0.)  #  [cm/day] "noFlux")#
    s.setBotBC("freeDrainage") #
    s.setParameter("Soil.BC.Top.SType", "2")  # michaelisMenten=8 (SType = Solute Type)
    s.setParameter("Soil.BC.Top.CValue", "0")  # michaelisMenten=8 (SType = Solute Type)
    s.setParameter("Soil.BC.Bot.SType", "6")  # michaelisMenten=8 (SType = Solute Type)
    s.setParameter("Soil.BC.Bot.CValue", "0.")
    s.setParameter( "Soil.css1Function", str(4))
    s.setParameter("Newton.EnableChop", "true")
    
    # UpwindWeight = 1, better when we have high solute gradient.
    # UpwindWeight = 0.5, better when have high water flow and low solute gradient
    s.setParameter("Flux.UpwindWeight", "1")#very important because we get high solute gradient.
    

    s.EnableResidualCriterion = False
    s.setParameter("Newton.EnableResidualCriterion", 
                     str( s.EnableResidualCriterion ))
    s.EnableAbsoluteResidualCriterion = False
    s.setParameter("Newton.EnableAbsoluteResidualCriterion", 
                     str( s.EnableAbsoluteResidualCriterion ))
    s.SatisfyResidualAndShiftCriterion = False
    s.setParameter("Newton.SatisfyResidualAndShiftCriterion",
                     str( s.SatisfyResidualAndShiftCriterion) )  
    s.MaxTimeStepDivisions = 10
    s.setParameter("Newton.MaxTimeStepDivisions",
                     str( s.MaxTimeStepDivisions) )  
    s.MaxSteps = 18
    s.setParameter("Newton.MaxSteps",
                     str( s.MaxSteps) )  
    s.MaxRelativeShift =1e-8
    s.setParameter("Newton.MaxRelativeShift", str(s.MaxRelativeShift))

    s.setParameter("Soil.IC.C", "0")  # g / cm3  # TODO specialised setter?

    s.setParameter("Component.MolarMass", "1.8e-2") 

    s.setParameter("Component.LiquidDiffusionCoefficient", "1.e-9")  # m2 s-1
    s.initializeProblem(maxDt = 0.01)

    # dummy sources of water and solutes to test the mass balance
    #source_map = { 0: -0.001, 1:-0.01, 2: 0.02, 3:-0.02, 4: 0.03, 5:-0.03}
    source_map = { 0: -0.01*0., 1:0.1, 2: 0.0, 3:0.0, 4: 0.0, 5:0.0}
    s.setSource(source_map)
    #source_map = { 0: 0.001, 1:-0.01, 2: 0.02, 3:-0.02, 4: 0.03, 5:-0.03}
    source_map = { 0: 0.001, 1:0.0, 2: 0.0, 3:0.0, 4: 0.0, 5:0.0}
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

        Wvolbefore = s.getCellVolumes() * s.getWaterContent() # cm3
        Smassbefore = s.getContent(1)#getSolution(1) * (Wvolbefore / 1.e3 ) # kg
        
        s.solve(dt)
        
        Wvolafter = s.getCellVolumes()*s.getWaterContent() # cm3
        Smassafter = s.getContent(1)# s.getSolution(1) * (Wvolafter / 1.e3 ) # kg
        
            
        scvFluxes = s.getFluxesPerCell(0) * dt  # cm3
        scvFluxesS = s.getFluxesPerCell(1) * dt  # kg
        scvSources = s.getSource(0) * s.getCellVolumes() * dt # cm3
        scvSourcesS = s.getSource(1) * s.getCellVolumes() * dt # kg
        
        if rank == 0:
            print('\tChange in water volume [cm3] per voxel:',Wvolafter-Wvolbefore)
            print('\tChange in solute mass [kg] per voxel:',Smassafter-Smassbefore)
            print('\n\n\tRMSE for water volume balance [cm3]:',np.mean(np.abs(scvFluxes+(Wvolafter-Wvolbefore)-scvSources)))
            print('\tRMSE for solute mass balance [kg]:',np.mean(np.abs(scvFluxesS+(Smassafter-Smassbefore)-scvSourcesS)),'\n\n')
            print(np.array(s.base.getScvfBoundaryFluxes()[0]))
            print(scvFluxes,scvSources)
        raise Exception



if __name__ == "__main__":


    simTimes = [0.5,0.78,1]  # days
    solve(simTimes)

