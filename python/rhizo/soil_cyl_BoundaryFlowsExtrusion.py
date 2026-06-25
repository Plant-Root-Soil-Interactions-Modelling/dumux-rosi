import sys; 
sys.path.append("../../experimental/fixedPointIter2/modules"); 
sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

from rosi_richards_cyl import RichardsCylFoam #as RichardsCylFoamInit # C++ part (Dumux binding)
from rosi_richards_cyl import RichardsCylFoamAna
#from rosi_richards_cylExtrusion import  RichardsCylFoamExtrusion # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part
from richards_no_mpi import RichardsNoMPIWrapper  

import matplotlib.pyplot as plt
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()
import numpy as np

""" 
Get the inner and outer boundary flows for the 1d axisymmetric soil
"""

r_in = 2  # cm
r_out = 10
length = 3

def solve(simtimes, N, analytical = False):

    loam = [0.08, 0.43, 0.04, 1.6, 5]  # K = 5 !
    if analytical:
        s = RichardsNoMPIWrapper(RichardsCylFoamAna()) 
    else:
        s = RichardsNoMPIWrapper(RichardsCylFoam()) 
        
    s.initialize()
    s.createGrid1d(np.linspace(r_in, r_out, N))  # [cm]
    s.setVGParameters([loam])

    s.setHomogeneousIC(-1000.)  # cm pressure head
    if analytical:
        s.setTopBC("constantFlux",0.1)  #  [cm/day] "noFlux")#
        s.setBotBC("constantFlux", -0.2) # "noFlux")#
        s.setParameter("Problem.useExtrusion", "true")
    else:
        s.setTopBC("constantFluxCyl",0.1)  #  [cm/day] "noFlux")#
        s.setBotBC("constantFluxCyl", -0.2) # "noFlux")#
        s.setParameter("Problem.useExtrusion", "false")
        
    s.initializeProblem(maxDt = 0.01)
    cellVolumes = s.getCellSurfacesCyl() * length # cm3
    if rank == 0:
        print(s)

    s.ddt = 1.e-3  # days

    simtimes.insert(0, 0)
    dt_ = np.diff(simtimes)
    
    watpots = []
    #x = s.get
    
    for r, dt in enumerate(dt_):
        

        time = simtimes[r] 
        print('time',time, s.useMoles())
            
        #if time >= 5:
        #    s.setSoluteTopBC([1], [0.])

        if rank == 0:
            print("*****", "#", r, "external time step", dt, " d, simulation time", s.simTime, "d, internal time step", s.ddt, "d")

        Wvolbefore = cellVolumes * s.getWaterContent() # cm3
        #Smassbefore = s.getSolution(1) * Wvolbefore # g
        print("Wvolbefore",Wvolbefore,sum(Wvolbefore))
        s.solve(dt, saveInnerDumuxValues_ = True)
        
        Wvolafter = cellVolumes*s.getWaterContent() # cm3
        print("Wvolafter",Wvolafter,sum(Wvolafter))
        #Smassafter = s.getSolution(1) * Wvolafter # g
        
            
        rootSoilFluxes = s.getInnerFlow(0, length) * dt # cm3
        print("rootSoilFluxes",rootSoilFluxes, s.getInnerFlow(0, length), dt)
        #rootSoilFluxesS = s.getInnerFlow(1, length) * dt # g
        soilSoilFluxes = s.getOuterFlow(0, length) * dt # cm3
        print("soilSoilFluxes",soilSoilFluxes, s.getOuterFlow(0, length), dt)
        print('getBoundaryFluxesPerFace_',s.getBoundaryFluxesPerFace_(0, length))
        print('getFaceSurfaces_',s.getFaceSurfaces_(length))
        #soilSoilFluxesS = s.getOuterFlow(1, length) * dt # g
        
        # TODO: currently, setSource not properly implemented for richards and richards 2c
        # so left out of the mass balance.
        # scvSources = s.getSource(0) * cellVolumes * dt # cm3
        # scvSourcesS = s.getSource(1) * cellVolumes * dt # kg
        
        if rank == 0:
            print('\tChange in water volume [cm3] per voxel:',sum(Wvolafter-Wvolbefore))
            #print('\tChange in solute mass [g] per voxel:',sum(Smassafter-Smassbefore))
            print('\tRMSE for water volume balance [cm3]:',np.mean(np.abs(rootSoilFluxes+soilSoilFluxes+sum(Wvolafter-Wvolbefore))))
            #print('\tRMSE for solute mass balance [g]:',np.mean(np.abs(rootSoilFluxesS+soilSoilFluxesS+sum(Smassafter-Smassbefore))),'\n\n')
            

       


if __name__ == "__main__":


    simTimes = [0.5,0.78,1]  # days
    simdataAna = solve(simTimes, 2, True)
    #simdata = solve(simTimes, 2, False)

"""
when we have no BC flow, seems to work ok, just inner flow.

todo:
- solve issue with bc flow, OK
- test source
- test analytical, OK
- weirdly sometime we do have flux at the outer boundary even though we said no flux
- att: i changed the dumux/assembly/cclocalresidual
"""