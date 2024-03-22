import sys; sys.path.append("./modules_fpit");  sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

from rosi_richardsnc import RichardsNCSP # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part

import matplotlib.pyplot as plt
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()
import numpy as np
import scenario_setup as scenario

""" 
Steady state linear solute transport in homogeneous soil profile
"""
"""
questions:
-y maxDt not respected
-y values different between nc and 10c <= implement constant vs millintonquirk diffusion
-at some point compaire again dumux and comsol for 10c
- y number of threads not always respected.
- att creation still uses all the threads
- after, still has control on all the trheads but only uses DUMUX_NUM_THREADS threads
"""

def solve(simtimes):

    loam = [0.08, 0.43, 0.04, 1.6, 5]  # K = 5 !
    s = RichardsWrapper(RichardsNCSP(),False)
    
    s.setParameter("Component.LiquidDiffusionCoefficient", "1.e-7")  # m2 s-1
    s.initialize()
    s.createGrid([-5., -5., -40.], [5., 5., 0.], [3,12,40])  # [cm]
    s.setVGParameters([loam])

    s.setHomogeneousIC(-5.)  # cm pressure head
    s.setTopBC("constantFlux", 0)  #  [cm/day]
    s.setBotBC("constantFlux", 0)

    s.setParameter("Soil.IC.C", "0")  # g / cm3  # TODO specialised setter?
    
    # TODO no idea, where this is neeeded, i don't want to use moles ever
    molarMass=12.011 #g/mol
    s.setParameter("Component.MolarMass", str(molarMass/1000)) #kg/mol  

    
#   
    soureC1 = 1e-4 *0#g/d
    s.setParameter("Soil.BC.Top.SType", "2") 
    s.setParameter("Soil.BC.Top.CValue", str(soureC1)) 
    s.setParameter("Soil.BC.Bot.SType", "2")  
    s.setParameter("Soil.BC.Bot.CValue", "0.")

    s.initializeProblem(maxDt = 250/(3600*24))

    s.wilting_point = -15000
    s.setCriticalPressure(s.wilting_point)  # for boundary conditions constantFlow, 
    

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
            s.ddt = 1.e-3  # days

            s.solve(1)
            
            if rank == 0:
                s.base.printParams()
            raise Exception

        x.append(s.getSolutionHead())
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

