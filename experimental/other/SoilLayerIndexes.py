import sys; sys.path.append("../../python/modules");  sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

from rosi_richardsnc import RichardsNCSP # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part

import matplotlib.pyplot as plt
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()

""" 
Steady state infiltration with a 3D SPGrid but no resolution in x and y (for speed), 
approximated with simulation time of one year

everything scripted, no input file needed

also works parallel with mpiexec
"""


def solve():
    soil = [0.045, 0.43, 0.15, 3, 1000]
    s = RichardsWrapper(RichardsNCSP())
    s.initialize()
    s.createGrid([-2., -2., -20.], [2., 2., 0.], [2, 2, 20])  # [cm]
    s.setHomogeneousIC(-50.)  # cm pressure head
    s.setBotBC("freeDrainage")
    s.setTopBC("constantFlux", 0)  #  [cm/day]
    s.setVGParameters([soil])
    s.initializeProblem()
    s.addVanGenuchtenDomain([-2., -2., -5.], [2., 2., 0.], 1)
    s.addVanGenuchtenDomain([-2., -2., -10.], [2., 2., -8.], 2)
    cellCenters = s.getCellCenters()
    inds = []
    for i in range(3):
        inds.append(s.getLayerCellIndx(i))
    if rank == 0:
        for i in range(3):
            print(f"index of cells in soil layer {i}:",inds[i])
            if len(inds[i]) > 0:
                print(f"centers of cells in soil layer {i}:",cellCenters[inds[i]])




if __name__ == "__main__":

    solve()

