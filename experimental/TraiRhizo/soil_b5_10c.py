import sys; sys.path.append("./modules_fpit"); 
sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

# RichardsNCSPnum
# RichardsNCSPSeq
# RichardsNCSPSSORC
# RichardsNCSPILURes
from rosi_richards10c import RichardsNCSPILU as RichardsNCSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part

import matplotlib.pyplot as plt
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()
import numpy as np
import scenario_setup as scenario
from cyl3plant_simple import simulate_const
""" 
Steady state linear solute transport in homogeneous soil profile
"""


def solve(simtimes):

    results_dir="./results/"

    # theta_r, theta_s, alpha, n, Ks
    loam = [0.08, 0.43, 0.04, 1.6, 5]  # K = 5 !
    usemoles = True
    s = RichardsWrapper(RichardsNCSP(), usemoles)
    s.soil =loam
    s.solidDensity = 2650 # [kg/m^3 solid] #taken from google docs TraiRhizo
    s.solidMolarMass = 60.08e-3 # [kg/mol] 
    s.setParameter( "Soil.MolarMass", str(s.solidMolarMass))
    s.setParameter( "Soil.solidDensity", str(s.solidDensity))
    s.css1Function=9
    s.setParameter( "Soil.css1Function", str(s.css1Function))
    noAds=True
    
    #scenario.setBiochemParam(s)
    s.setParameter("1.Component.LiquidDiffusionCoefficient", "1.e-7")
    #s.setParameter("Assembly.Multithreading","false")
    s.initialize()#doMPI_=False)
    s.createGrid([-5., -5., -40.], [5., 5., 0.], [3,12,40])  # [cm]
    s.setVGParameters([loam])

    s.setHomogeneousIC(-5.)  # cm pressure head
    s.setTopBC("constantFlux", 0)  #  [cm/day]
    s.setBotBC("constantFlux", 0)

    for i in range(s.numComp):
        s.setParameter( "Soil.IC.C"+str(i+1), str(0 ))
        
    
    
    molarMass=12.011 #g/mol
    soureC1 = 1e-4 #g/d
    s.setParameter( "Soil.BC.Top.C1Type", str(2)) 
    s.setParameter( "Soil.BC.Top.C1Value", str(0))#soureC1/molarMass)) 
    s.setParameter( "Soil.BC.Bot.C1Type", str(2)) 
    s.setParameter( "Soil.BC.Bot.C1Value", str(0)) 
    #s.initializeProblem(maxDt = 250/(3600*24))
    for i in range(2, s.numComp+1):
        s.setParameter( "Soil.BC.Bot.C"+str(i)+"Type", str(2))
        s.setParameter( "Soil.BC.Top.C"+str(i)+"Type", str(2))
        s.setParameter( "Soil.BC.Bot.C"+str(i)+"Value", str(0)) 
        s.setParameter( "Soil.BC.Top.C"+str(i)+"Value", str(0 )) 
        
        
    
        
    s.initializeProblem(maxDt = 250/(3600*24))
    
        
    s.wilting_point = -15000
    s.setCriticalPressure(s.wilting_point)  # for boundary conditions 
    
    x, c = [], []
    simtimes.insert(0, 0)
    dt_ = np.diff(simtimes)
    s.ddt = 0.1
    for r, dt in enumerate(dt_):

        for i in range(0, int(dt)):

            time = simtimes[r] + i
            print('time',time)

            #if time >= 5:
            #    s.setSoluteTopBC([1], [0.])

            if rank == 0:
                print("*****", "#", i, "external time step", dt, " d, simulation time", s.simTime, "d, internal time step", s.ddt, "d")

            #simulate_const(s,1)
            s.solve(1)#, doMPIsolve_=True)
                    
            #if rank == 0:
            #    s.base.printParams()
        x.append(s.getSolutionHead())
        #print(max(s.getSolution(1)))
        c.append(s.getSolution(1))

    points = s.getDofCoordinates()

    return x, c, points


if __name__ == "__main__":

    cols = ["r*", "g*", "b*", "m*"] * 10
    simtimes = [5, 10, 20, 30]  # days
    xa, ca, za = solve([5, 10, 20, 30])
    # xb, cb, zb, = solve([1.])
    # xc, cc, zc = solve([10.])

    if rank == 0 and False:
        # plt.plot(xa, za[:, 2], "r*")
        for i in range(0, len(simtimes)):
            plt.plot(ca[i], za[:, 2], cols[i])
        plt.legend(simtimes)

        # plt.plot(RichardsWrapper.to_head(xb), zb[:, 2], "g*")
        # plt.plot(RichardsWrapper.to_head(xc), zc[:, 2], "b*")
        plt.show()
        for i in range(0, len(simtimes)):
            plt.plot(xa[i], za[:, 2], cols[i])
        plt.show()

