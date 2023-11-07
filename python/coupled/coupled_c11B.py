import sys; sys.path.append("../modules_fpit/"); sys.path.append("../../../CPlantBox/");  sys.path.append("../../../CPlantBox/src/")
sys.path.append("../../build-cmake/cpp/python_binding/")

#import visualisation.vtk_plot as vp
from functional.xylem_flux import XylemFluxPython
#from xylem_flux import XylemFluxPython
import plantbox as pb
import rsml.rsml_reader as rsml
from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part

from math import *
import numpy as np
import matplotlib.pyplot as plt
import timeit
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()

""" 
Benchmark M1.1 single root in in soil (using the classic sink)

additionally, compares exact root system flux (Meunier et al.) with approximation

also works parallel with mpiexec (only slightly faster, due to overhead)
"""

r_root = 0.02  # cm
r_out = 0.6

        
        
def write_file_array(name, data, space =",", directory_ ="./results/", fileType = '.csv', allranks = False ):
    print('write_file_array', name, directory_)
    if (rank == 0) or allranks:
        name2 = directory_+ name+ fileType
        with open(name2, 'a') as log:
            log.write(space.join([num for num in map(str, data)])  +'\n')

def solve(soil, simtimes, q_r, N):
    """ 
    soil            soil type
    simtimes        simulation output times
    q_r             root flux [cm/day]   
    N               spatial resolution
    """

    l = 0.53#np.sqrt((r_out * r_out - r_root * r_root) * np.pi) / 2  # same area as cylindrical

    """ Soil problem """
    s = RichardsWrapper(RichardsSP(), False)
    s.initialize()
    s.setHomogeneousIC(-100)  # cm pressure head
    s.setTopBC("noflux")
    s.setBotBC("noflux")
    print('to create grid')
    s.createGrid([-l, -l, -1.], [l, l, 0.], [N, N, 1])  # [cm]
    print('grid created')
    s.setVGParameters([soil[0:5]])
    print('to initializeProblem')
    s.initializeProblem()
    print('initializeProblem created')
    s.setCriticalPressure(-15000)
    print('getGlobal2localCellIdx', s.base.getGlobal2localCellIdx())
    print('getLocal2globalPointIdx', s.base.getLocal2globalPointIdx())

    return 0

if __name__ == "__main__":

    sand = [0.045, 0.43, 0.15, 3, 1000, "Sand"]
    loam = [0.08, 0.43, 0.04, 1.6, 50, "Loam"]
    clay = [0.1, 0.4, 0.01, 1.1, 10, "Clay"]

    sim_times = np.linspace(0, 25, 250)  # temporal resolution of 0.1 d

    if rank == 0:
        fig, ax = plt.subplots(2, 3, figsize=(14, 14))
        t0 = timeit.default_timer()

    jobs = [loam, 0.1, 0, 1]#[sand, 0.1, 0, 0], [loam, 0.1, 0, 1])#, [clay, 0.1, 0, 2], [sand, 0.05, 1, 0], [loam, 0.05, 1, 1], [clay, 0.05, 1, 2])
    soil, qj, i, j = jobs
    print('to next soil', soil, qj, i, j)
    y, x, t = solve(soil, sim_times, qj, 10)
    print('did soil')
    #print(jobs)
    #soil, qj, i, j = jobs
    #print(soil, qj, i, j)
    #soil, qj, i, j = jobs
    #y, x, t = solve(soil, sim_times, qj, 41)
    


