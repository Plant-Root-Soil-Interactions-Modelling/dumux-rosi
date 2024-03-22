import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

from rosi_richardsnc import RichardsNCSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part

import matplotlib.pyplot as plt
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()
import numpy as np

""" 
Steady state linear solute transport in homogeneous soil profile
"""


def solve(simtimes):

    loam = [0.08, 0.43, 0.04, 1.6, 5]  # K = 5 !
    s = RichardsWrapper(RichardsNCSP())

    s.initialize()
    s.createGrid([-5., -5., -200.], [5., 5., 0.], [1, 1, 1990])  # [cm]
    s.setVGParameters([loam])

    # theta = 0.378, benchmark is set be nearly fully saturated, so we don't care too much about the specific values
    s.setHomogeneousIC(-5.)  # cm pressure head
    s.setTopBC("constantFlux", 2)  #  [cm/day]
    s.setBotBC("freeDrainage")

    s.setParameter("Soil.IC.C", "0")  # g / cm3  # TODO specialised setter?

    s.setParameter("Component.MolarMass", "1.8e-2")  # TODO no idea, where this is neeeded, i don't want to use moles ever

    s.setParameter("Component.LiquidDiffusionCoefficient", "1.e-9")  # m2 s-1
#     s.setParameter("SpatialParams.Tortuosity", "0.001")  # porosity * 0.3

    # s.setSTopBC([1], [1.e-4])  # 1 == Dirchlet, change for t>1 days
    # s.setSBotBC([6], [0.])  # 6 == outflow
    s.setParameter("Soil.BC.Top.SType", "2")  # michaelisMenten=8 (SType = Solute Type)
    s.setParameter("Soil.BC.Top.CValue", "1.e-4")  # michaelisMenten=8 (SType = Solute Type)
    s.setParameter("Soil.BC.Bot.SType", "6")  # michaelisMenten=8 (SType = Solute Type)
    s.setParameter("Soil.BC.Bot.CValue", "0.")

    s.initializeProblem(maxDt = 0.01)

    if rank == 0:
        print(s)

    x, c = [], []
    s.ddt = 1.e-3  # days

    simtimes.insert(0, 0)
    dt_ = np.diff(simtimes)
    for r, dt in enumerate(dt_):

        for i in range(0, int(dt)):

            time = simtimes[r] + i
            print('time',time)

            #if time >= 5:
            #    s.setSoluteTopBC([1], [0.])

            if rank == 0:
                print("*****", "#", i, "external time step", dt, " d, simulation time", s.simTime, "d, internal time step", s.ddt, "d")

            s.solve(1)

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

    if rank == 0:
        # plt.plot(xa, za[:, 2], "r*")
        for i in range(0, len(simtimes)):
            plt.plot(ca[i], za[:, 2], cols[i])
        plt.legend(simtimes)

        # plt.plot(RichardsWrapper.to_head(xb), zb[:, 2], "g*")
        # plt.plot(RichardsWrapper.to_head(xc), zc[:, 2], "b*")
        plt.show()

