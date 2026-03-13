import sys; sys.path.append("../../experimental/fixedPointIter2/modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

# from rosi_richardsnc_cyl import RichardsNCCylFoam  # C++ part (Dumux binding)
# from richards import RichardsWrapper  # Python part

from rosi_richards4c_cyl import Richards4CCylFoam as RichardsNCCylFoam# C++ part (Dumux binding)
from richards_no_mpi import RichardsNoMPIWrapper  # Python part of cylindrcial model (a single cylindrical model is not allowed to run in parallel)
from rosi_richards4c import Richards4CSPILU as Richards4CSP  # C++ part (Dumux binding), macroscopic soil model
from richards import RichardsWrapper  # Python part, macroscopic soil model

import matplotlib.pyplot as plt
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()
import numpy as np

""" 
Get the inner and outer boundary flows for the 1d axisymmetric soil
"""

r_in = 2  # cm
r_out = 10
length = 3


solidDensity = 2650 # [kg/m^3 solid] #taken from google docs TraiRhizo
solidMolarMass = 60.08e-3 # [kg/mol] 
theta_S = 0.494
# [mol / m3 solid] =[kg/m^3 solid] / [kg/mol] 
solidMolDensity = solidDensity/solidMolarMass
# [mol / m3 scv] = [mol / m3 solid] * [m3 solid /m3 space]
bulkDensity_m3 = solidMolDensity*(1.- theta_S)
bulkMassDensity_gpercm3 = solidDensity*(1.- theta_S)*1000/1e6


soil_type = 0
kads = 7.07e+02 # m3/kgC/yr, see 10.1016/j.soilbio.2020.107912, A.3
kdes =  1.63e+03 # [1/yr] see 10.1016/j.soilbio.2020.107912, A.3
k_clay_silt = {}
k_clay_silt[0] = 0.67
k_clay_silt[1] = 0.082
molarMassC = 12.011
mg_per_molC = molarMassC * 1000.
yr_per_d = 1/365 # [yr/d]
m3_per_cm3 = 1e-6; # m3/cm3
cm3_per_m3 = 1e6; # cm3/m3

# [kg/g] * [g/mol] = kg/mol
kgC_per_mol = (1/1000) * molarMassC
# [m3/kgC/yr] * [yr/d] * [cm3/m3] * [kgC/mol] = [cm3/mol/d]
kads = kads * yr_per_d * cm3_per_m3 * kgC_per_mol # [cm3/mol/d]
kdes = kdes * yr_per_d # [1/d]
Qmmax = k_clay_silt[soil_type] * 0.079 # max ratio gOC-gmineral soil, see 10.1016/j.soilbio.2020.107912
# [g OC / g mineral soil] * [g mineral soil/ cm3 bulk soil] *[ mol C/g C]
CSSmax_ = Qmmax * bulkMassDensity_gpercm3*(1/molarMassC)
CSSmax = CSSmax_ # mol C/cm3 bulk soil


def solve(simtimes, N):

    loam = [0.08, 0.43, 0.04, 1.6, 5]  # K = 5 !
    s = RichardsWrapper(RichardsNCCylFoam())

    s.initialize()
    lb = 0.5
    points = np.logspace(np.log(r_in) / np.log(lb), np.log(r_out) / np.log(lb), 
                         10, base = lb)

    s.createGrid1d(points)# cm
                
    s.setVGParameters([loam])

    # theta = 0.378, benchmark is set be nearly fully saturated, so we don't care too much about the specific values
    s.setHomogeneousIC(-1000.)  # cm pressure head

    s.setParameter("Problem.segLength", str(1))  # cm
    s.setParameter("SpatialParams.Temperature","293.15") # todo: redefine at each time step
    s.setParameter("Soil.BC.dzScaling", "1")
    s.setParameter( "Soil.css1Function", str(9))
    s.setParameter("Problem.verbose", "0")
    s.setParameter("Problem.reactionExclusive", "0")
    s.setParameter("Soil.MolarMass", str(solidMolarMass))
    s.setParameter("Soil.solidDensity", str(solidDensity))
    s.setParameter("Flux.UpwindWeight", "1")


    s.setParameter("Soil.kads", str(1)) #[cm3/mol/d]
    s.setParameter("Soil.kdes", str(1)) #[1/d]            
    s.setParameter("Soil.CSSmax", str(0)) #[mol/cm3 scv zone 1] or mol
    s.setParameter("Soil.vmax_decay", str(0)) #mol C / m^3 scv / s
    s.setParameter("Problem.doDecay",str(0))

            
    s.setTopBC("constantFluxCyl",0.1)  #  [cm/day] "noFlux")#
    s.setBotBC("constantFluxCyl", -0.2) # "noFlux")#
    #s.setParameter("Soil.BC.Top.SType", "3")  # michaelisMenten=8 (SType = Solute Type)
    #s.setParameter("Soil.BC.Top.CValue", "-0.05")  # michaelisMenten=8 (SType = Solute Type)
    #s.setParameter("Soil.BC.Bot.SType", "3")  # michaelisMenten=8 (SType = Solute Type)
    #s.setParameter("Soil.BC.Bot.CValue", "0.1")
    s.setParameter("Soil.IC.C1", cyl.dumux_str(1) ) 
    s.setParameter("Soil.IC.C", "0")  # g / cm3  # TODO specialised setter?

    s.setParameter("Component.MolarMass", "1.8e-2") 

    s.setParameter("Component.LiquidDiffusionCoefficient", "1.e-9")  # m2 s-1
    s.initializeProblem(maxDt = 0.01)

    #cellVolumes = s.getCellSurfacesCyl() * length # cm3
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

        Wvolbefore = s.getWaterVolumes() # cellVolumes * s.getWaterContent() # cm3
        Smassbefore = s.getSolution(1) * Wvolbefore # g
        
        s.solve(dt, saveInnerFluxes_ = True)
        
        Wvolafter =s.getWaterVolumes() # cellVolumes*s.getWaterContent() # cm3
        Smassafter = s.getSolution(1) * Wvolafter # g
        
            
        rootSoilFluxes = s.getInnerFlow(0, length) * dt # cm3
        rootSoilFluxesS = s.getInnerFlow(1, length) * dt # g
        soilSoilFluxes = s.getOuterFlow(0, length) * dt # cm3
        soilSoilFluxesS = s.getOuterFlow(1, length) * dt # g
        
        # TODO: currently, setSource not properly implemented for richards and richards 2c
        # so left out of the mass balance.
        # scvSources = s.getSource(0) * cellVolumes * dt # cm3
        # scvSourcesS = s.getSource(1) * cellVolumes * dt # kg
        
        if rank == 0:
            print('\tChange in water volume [cm3] per voxel:',sum(Wvolafter-Wvolbefore))
            print('\tChange in solute mass [g] per voxel:',sum(Smassafter-Smassbefore))
            print('\n\n\tRMSE for water volume balance [cm3]:',np.mean(np.abs(rootSoilFluxes+soilSoilFluxes+sum(Wvolafter-Wvolbefore))))
            print('\tRMSE for solute mass balance [g]:',np.mean(np.abs(rootSoilFluxesS+soilSoilFluxesS+sum(Smassafter-Smassbefore))),'\n\n')
            

       


if __name__ == "__main__":


    simTimes = [0.5,0.78,1]  # days
    solve(simTimes, 2)

