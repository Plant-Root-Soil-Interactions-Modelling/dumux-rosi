import sys; sys.path.append("../modules_fpit/"); sys.path.append("../../../CPlantBox/");  sys.path.append("../../../CPlantBox/src/")
sys.path.append("../../build-cmake/cpp/python_binding/")

from functional.xylem_flux import XylemFluxPython  # Python hybrid solver
import plantbox as pb
import rsml.rsml_reader as rsml
from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part

from math import *
import numpy as np
import matplotlib.pyplot as plt
import timeit
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size()

""" 
Benchmark M1.1 single root in in soil (using the classic sink)

additionally, compares exact root system flux (Meunier et al.) with approximation

also works parallel with mpiexec (only slightly faster, due to overhead)
"""

def write_file_float(name, data, directory_, allranks = False):
    if (rank == 0) or allranks:
        name2 = directory_+ name+ '.txt'
        with open(name2, 'a') as log:
            log.write(repr( data)  +'\n')
        
def write_file_array(name, data, space =",", directory_ ="./results/", fileType = '.txt', allranks = False ):
    if (rank == 0) or allranks:
        name2 = directory_+ name+ fileType
        with open(name2, 'a') as log:
            log.write(space.join([num for num in map(str, data)])  +'\n')
r_root = 0.02  # cm
r_out = 0.6


def solve(soil, simtimes, q_r, N):
    """ 
    soil            soil type
    simtimes        simulation output times
    q_r             root flux [cm/day]   
    N               spatial resolution
    """

    q_r = q_r * 1  # * 1 g/cm3 = g/cm2/day
    q_r = q_r * 2 * r_root * np.pi * 1.  # g/day
    print("Qr as sink", q_r, "g day-1")
    trans = q_r
    l = np.sqrt((r_out * r_out - r_root * r_root) * np.pi) / 2  # same area as cylindrical

    """ Root problem"""
    n1 = pb.Vector3d(0, 0, 0)
    n2 = pb.Vector3d(0, 0, -0.5)
    n3 = pb.Vector3d(0, 0, -1)
    seg1 = pb.Vector2i(0, 1)
    seg2 = pb.Vector2i(1, 2)
    #  rs = pb.MappedSegments([n1, n3], [seg1], [r_root])  # a single root
    rs = pb.MappedSegments([n1, n2, n3], [seg1, seg2], [r_root, r_root])  # a single root with 2 segments
    rs.setRectangularGrid(pb.Vector3d(-l, -l, -1.), pb.Vector3d(l, l, 0.), pb.Vector3d(N, N, 1), False)
    r = XylemFluxPython(rs)
    r.setKr([1.e-7])
    r.setKx([1.e-7])

    """ Soil problem """
    s = RichardsWrapper(RichardsSP(), False)
    s.initialize()
    print('did init')
    s.setHomogeneousIC(-100)  # cm pressure head
    print('w1')
    s.createGrid([-l, -l, -1.], [l, l, 0.], [N, N, 1])  # [cm]
    print('w2')
    s.setTopBC("noflux")
    print('w3')
    s.setBotBC("noflux")
    print('w4')
    s.setVGParameters([soil[0:5]])
    print('w5')
    s.initializeProblem()
    print('w6')
    s.setCriticalPressure(-15000)
    print('crit pressure')
    """ Coupling (map indices) """
    picker = lambda x, y, z: s.pick([x, y, z])
    r.rs.setSoilGrid(picker)
    cci = picker(0, 0, 0)  # collar cell index
    print('picker')
    """ Numerical solution """
    start_time = timeit.default_timer()
    nsp = N  # number of sample points for output
    y_ = np.linspace(0, l, nsp)
    y_ = np.expand_dims(y_, 1)
    x_ = np.hstack((np.zeros((nsp, 1)), y_, np.zeros((nsp, 1))))

    dt_ = np.diff(simtimes)
    s.ddt = 1.e-5  # initial Dumux time step [days]
    print('to dt_')
    for dt in dt_:
        sx = s.getSolutionHead()
        print('sx')
        if rank == 0:  # Root part is not parallel
            print('solve_neumann, cci', cci, 'trans', trans)
            rx = r.solve_neumann(0., -trans, sx, True)  # xylem_flux.py
            print('solve_neumann finished',min(rx), max(rx) )
            fluxes_approx = r.soilFluxes(0., rx, sx, True)  # class XylemFlux is defined in MappedOrganism.h
            print('fluxes_approx',fluxes_approx)
            fluxes_exact = r.soilFluxes(0., rx, sx, False)  # class XylemFlux is defined in MappedOrganism.h
            print('fluxes_exact',fluxes_exact)
            collar_flux = r.collar_flux(0., rx, sx)
            print('collar_flux',collar_flux)
            off = abs(100 * (1 - fluxes_approx[cci] / (-trans)))
            print('off',off)
            print('cci, fluxes_approx[cci], fluxes_exact[cci], collar_flux[0], off',cci, fluxes_approx[cci], fluxes_exact[cci], collar_flux, off)
            print("fluxes at {:g}: approx {:g}, exact {:g}, collar flux {:g} [g day-1], approximation is {:g}% off"
                  .format(cci, fluxes_approx[cci], fluxes_exact[cci], collar_flux, off))
            fluxes = fluxes_exact
            sum_flux = 0.
            for f in fluxes.values():
                sum_flux += f
            print("Summed fluxes {:g} = {:g} [g day-1],  index 0: {:g}\n".format(sum_flux, -trans, fluxes[cci]))
        else:
            fluxes = None
            rx = None
        print('share data')
        fluxes = comm.bcast(fluxes, root = 0)  # Soil part runs parallel
        rx = comm.bcast(rx, root = 0)  # Soil part runs parallel
        print('set source')
        s.setSource(fluxes)  # g day-1, richards.py
        print('did set source')
        write_file_array('rx'+str(max_rank),rx, directory_ ="./", fileType = '.csv')
        write_file_float('fluxes'+str(max_rank),fluxes, directory_ ="./")
        print("to solve")
        s.solve(dt)
        print("finished solve")
        x0 = s.getSolutionHeadAt(cci)
        if x0 < -15000:
            if rank == 0:
                print("Simulation time at -15000 cm > {:.3f} cm after {:.3f} days".format(float(x0), s.simTime))
            testrx = np.array(s.interpolateNN(x_))
            write_file_array('testrx'+str(max_rank),testrx, directory_ ="./", fileType = '.csv')
            print('np.array(s.interpolateNN(x_))',rank,np.array(s.interpolateNN(x_)))
            y = s.to_head(np.array(s.interpolateNN(x_)))
            return y, x_[:, 1], s.simTime

    y = s.toHead(s.interpolateNN(x_))
    return y, x_[:, 1], s.simTime


if __name__ == "__main__":

    sand = [0.045, 0.43, 0.15, 3, 1000, "Sand"]
    loam = [0.08, 0.43, 0.04, 1.6, 50, "Loam"]
    clay = [0.1, 0.4, 0.01, 1.1, 10, "Clay"]

    sim_times = np.linspace(0, 25, 250)  # temporal resolution of 0.1 d

    jobs = ( [loam, 0.1, 0, 1], [clay, 0.1, 0, 2], [sand, 0.05, 1, 0], [loam, 0.05, 1, 1], [clay, 0.05, 1, 2])
    for soil, qj, i, j in jobs:
        y, x, t = solve(soil, sim_times, qj, 41)
        write_file_array('soilHeadnew'+str(max_rank),y, directory_ ="./", fileType = '.csv')

    if rank == 0:
        print("Elapsed time: ", timeit.default_timer() - t0, "s")
        plt.show()

