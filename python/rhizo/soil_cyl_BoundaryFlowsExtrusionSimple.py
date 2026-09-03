import sys; 
sys.path.append("../modules"); 
sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

from rosi_richards_cyl import RichardsCylFoam #as RichardsCylFoamInit # C++ part (Dumux binding)
from rosi_richards_cyl import RichardsCylFoamAna
#from rosi_richards_cylExtrusion import  RichardsCylFoamExtrusion # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part
from richards_no_mpi_flat import RichardsNoMPIFlatWrapper 
import functional.van_genuchten as vg

import matplotlib.pyplot as plt
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()
import numpy as np

""" 
Get the inner and outer boundary flows for the 1d axisymmetric soil
"""

r_in = 0.01  # cm
r_out = 1
length = 1
    
def solve(simtimes, N, analytical = False):

    #loam = [0.08, 0.43, 0.04, 1.6, 5]  # K = 5 !
    if analytical:
        s = RichardsNoMPIFlatWrapper(RichardsCylFoamAna()) 
    else:
        s = RichardsNoMPIFlatWrapper(RichardsCylFoam()) 
        
        
    #s.setParameter( "Assembly.NumericDifference.BaseEpsilon", str(1e-10));
    #s.setParameter("Assembly.NumericDifference.PriVarMagnitude", s.dumux_str([1e-10, 1e-2, 1e-10, 1e-10]));
    s.initialize()
    s = setDefault(s)
    s.doAds = False
    s.doDecay = False
    s.doBioChemicalReaction = False
    s.setParameter("Problem.segLength", str(length))#
    s.setParameter("Soil.MolarMass", str(length))#
    s.createGrid1d(np.linspace(r_in, r_out, N))  # [cm]
    #s.setVGParameters([loam])
    soilTextureAndShape = getSoilTextureAndShape()
    setSoilParam(s, soilTextureAndShape)
    getBiochemParam(s,soil_type = 0,sorp = 0, diff = 0)
    setBiochemParam(s) 

    s.setHomogeneousIC(-1000.)  # cm pressure head
    s.setTopBC("constantFlux", 0.1)  #  [cm/day] "noFlux")#
    s.setBotBC("constantFlux", -0.2) # "noFlux")#
    cin = 9.8703454e-07 # mol /day
    for i in range(1, 2):#s.numComp):# no flux
        s.setParameter( "Soil.BC.Bot.C"+str(i)+"Type", str(2))
        s.setParameter( "Soil.BC.Bot.C"+str(i)+"Value", str(cin / (2 * np.pi * r_in * length))) # mol / cm -2 /day
    s.setParameter( "Soil.IC.C1", str(0 ))
        
    if analytical:
        s.setParameter("Problem.useExtrusion", "true")
    else:
        s.setParameter("Problem.useExtrusion", "false")
        
    s.initializeProblem(maxDt = 0.01)
    #s.eps_regularization = 1e-10 
    #s.setRegularisation(s.eps_regularization, s.eps_regularization) 
    cellVolumes = s.getCellSurfacesCyl() * length # cm3
    if rank == 0:
        print(s)

    s.ddt = 1.e-3  # days

    simtimes.insert(0, 0)
    dt_ = np.diff(simtimes)
    
    watpots = []
    #x = s.get
    withsol = False
    for r, dt in enumerate(dt_):
        

        time = simtimes[r] 
        print('time',time, s.useMoles)
            
        #if time >= 5:
        #    s.setSoluteTopBC([1], [0.])

        if rank == 0:
            print("*****", "#", r, "external time step", dt, " d, simulation time", s.simTime, "d, internal time step", s.ddt, "d")

        Wvolbefore = cellVolumes * s.getWaterContent() # cm3
        print("Wvolbefore",Wvolbefore,sum(Wvolbefore))
        if withsol:
            Smassbefore = s.getContent(1) # mol  * Wvolbefore
            print("Smassbefore",Smassbefore,sum(Smassbefore))
        s.solve(dt, saveInnerFluxes_ = True)
        
        Wvolafter = cellVolumes*s.getWaterContent() # cm3
        print("Wvolafter",Wvolafter,sum(Wvolafter))
        if withsol:
            Smassafter = s.getContent(1) #* Wvolafter # mol
            print("Smassafter",Smassafter,sum(Smassafter))
        
            
        rootSoilFluxes = s.getInnerFlow(0, length) * dt # cm3
        if withsol:
            rootSoilFluxesS = s.getInnerFlow(1, length) * dt # mol
        soilSoilFluxes = s.getOuterFlow(0, length) * dt # cm3
        if withsol:
            soilSoilFluxesS = s.getOuterFlow(1, length) * dt # mol
        print("rootSoilFluxes",rootSoilFluxes, s.getInnerFlow(0, length), dt)
        print("soilSoilFluxes",soilSoilFluxes, s.getOuterFlow(0, length), dt)
        print('getBoundaryFluxesPerFace_',s.getBoundaryFluxesPerFace_(0, length)) # mol/day
        if withsol:
            print('getFaceSurfaces_',s.getFaceSurfaces_(length)) # cm2
            print("base.getScvfBoundaryFluxes()",s.base.getScvfBoundaryFluxes()[1][0] *24 * 3600 / 1e4) # mol/cm2/day
            print('neumann', s.base.getNeumann(0,1))
        
        # TODO: currently, setSource not properly implemented for richards and richards 2c
        # so left out of the mass balance.
        # scvSources = s.getSource(0) * cellVolumes * dt # cm3
        # scvSourcesS = s.getSource(1) * cellVolumes * dt # kg
        
        if rank == 0:
            print('\tChange in water volume [cm3] per voxel:',sum(Wvolafter-Wvolbefore))
            if withsol:
                print('\tChange in solute mass [g] per voxel:',sum(Smassafter-Smassbefore))
            print('\tRMSE for water volume balance [cm3]:',np.mean(np.abs(rootSoilFluxes+soilSoilFluxes+sum(Wvolafter-Wvolbefore))))
            if withsol:
                print('\tRMSE for solute mass balance [g]:',np.mean(np.abs(rootSoilFluxesS+soilSoilFluxesS+sum(Smassafter-Smassbefore))),'\n\n')
        #print('getContent',s.getContent(1), s.getContent(2), s.getContent(3))   

       


if __name__ == "__main__":


    simTimes = [0.5]#,0.78,1]  # days
    #simdataAna = solve(simTimes, 3, True)
    simdata = solve(simTimes, 3, False)

"""
when we have no BC flow, seems to work ok, just inner flow.

todo:
- Y results so different for solutes in ana and numeric?
- error when i have more than one cell
- solutes works on its own but not when water is moving
    - jacobian not well filled?
    - cach not well rolled back or updated?
    - ussie that carbon advection according to numercial diff vs water flow from analytical diff?
    - ok so now inner BC is wrong for the numerical
"""