import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

from functional.xylem_flux import XylemFluxPython  # Python hybrid solver
import plantbox as pb
import rsml.rsml_reader as rsml
from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from rosi_richards_cyl import RichardsCylFoam  # C++ part (Dumux binding)
from richards import RichardsWrapper  # for soil part wrappring RichardsSP
from richards_no_mpi import RichardsNoMPIWrapper  # for cylindrical part, RichardsCylFoam must be compiled disabeling MPI

from math import *
import numpy as np
import matplotlib.pyplot as plt
import timeit
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()

""" 
Benchmark M1.1 single root in in soil

Prestudy for cylindrical models around segments:
soil sink is based on the cylindrical model, no information from mascroscopic soil to cylindrical model (yet)

also works parallel with mpiexec (only slightly faster, due to overhead)
"""

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
    rs = pb.MappedSegments([n1, n3], [seg1], [r_root])  # a single root
    # rs = pb.MappedSegments([n1, n2, n3], [seg1, seg2], [r_root, r_root])  # a single root with 2 segments
    rs.setRectangularGrid(pb.Vector3d(-l, -l, -1.), pb.Vector3d(l, l, 0.), pb.Vector3d(N, N, 1), False)
    r = XylemFluxPython(rs)
    r.setKr([1.e-7])
    r.setKx([1.e-7])

    """ Soil problem """
    s = RichardsWrapper(RichardsSP())
    s.initialize()
    s.setHomogeneousIC(-100)  # cm pressure head
    s.setTopBC("noflux")
    s.setBotBC("noflux")
    s.createGrid([-l, -l, -1.], [l, l, 0.], [N, N, 1])  # [cm]
    s.setVGParameters([soil[0:5]])
    s.initializeProblem()
    s.setCriticalPressure(-15000)
    # s.setRegularisation(1.e-6, 1.e-6)

    """ Cylindrical problems """
    rich_cyls = []
    ns = len(r.rs.segments)
    rsx = np.zeros((ns,))  # root soil interface
    for i in range(0, ns):
        rcyl = RichardsNoMPIWrapper(RichardsCylFoam())
        rcyl.initialize([""], False)  # [""], False
        rcyl.createGrid([r_root], [r_out], [N])  # [cm]
        rcyl.setHomogeneousIC(-100.)  # cm pressure head
        rcyl.setVGParameters([soil[0:5]])
        rcyl.setOuterBC("fluxCyl", 0.)  #  [cm/day]
        rcyl.setInnerBC("fluxCyl", 0.)  # [cm/day]
        rcyl.initializeProblem()
        rcyl.setCriticalPressure(-15000)  # cm pressure head
        rich_cyls.append(rcyl)

    """ Coupling (map indices) """
    picker = lambda x, y, z: s.pick([x, y, z])
    r.rs.setSoilGrid(picker)
    cci = picker(0, 0, 0)  # collar cell index

    """ Numerical solution """
    start_time = timeit.default_timer()
    nsp = N  # number of sample points for output
    y_ = np.linspace(0, l, nsp)
    y_ = np.expand_dims(y_, 1)
    x_ = np.hstack((np.zeros((nsp, 1)), y_, np.zeros((nsp, 1))))

    dt_ = np.diff(simtimes)
    s.ddt = 1.e-5  # initial Dumux time step [days]

    for dt in dt_:

        sx = s.getSolutionHead()

        if rank == 0:  # Root part is not parallel

            for j, rc in enumerate(rich_cyls):
                rsx[j] = rc.getInnerHead()

            rx = r.solve_neumann(0., -trans, rsx, False)  # xylem_flux.py, True: soil matric potential per cell, False: soil matric potential per segment

            # For the soil model
            seg_fluxes = r.segFluxes(0., rx, rsx, approx = False)  # class XylemFlux is defined in CPlantBox XylemFlux.h
            soil_fluxes = r.sumSegFluxes(seg_fluxes)  # class XylemFlux is defined in CPlantBox XylemFlux.h

            fluxes_exact = r.soilFluxes(0., rx, sx, False)  # class XylemFlux is defined in MappedOrganism.h
            collar_flux = r.collar_flux(0., rx, sx)
            # print("fluxes at {:g}, exact {:g}, collar flux {:g} [g day-1]".format(cci, fluxes_exact[cci], collar_flux[0]))
            print("fluxes at {:g},  soil {:g}".format(cci, soil_fluxes[cci]))

            fluxes = soil_fluxes
            sum_flux = 0.
            for f in fluxes.values():
                sum_flux += f
            print("Summed fluxes {:g}".format(sum_flux + trans))

            # run cylindrical model
            for j, rc in enumerate(rich_cyls):  # set sources
                l = 1.
                print("Set inner flux to {:g} [cm day-1]\n".format(seg_fluxes[j] / (2 * np.pi * r_root * l)))
                rsx[j] = rc.setInnerFluxCyl(seg_fluxes[j] / (2 * np.pi * r_root * l))  # /
            for j, rc in enumerate(rich_cyls):  # simualte time step
                rc.solve(dt)
        else:
            fluxes = None

        fluxes = comm.bcast(fluxes, root = 0)  # Soil part runs parallel
        s.setSource(fluxes)  # g day-1, richards.py

        s.solve(dt)

        x0 = s.getSolutionHeadAt(cci)
        if x0 < -15000:
            if rank == 0:
                print("Simulation time at -15000 cm > {:.3f} cm after {:.3f} days".format(float(x0), s.simTime))
            y = s.to_head(np.array(s.interpolateNN(x_)))
            yc = rich_cyls[-1].getSolutionHead()
            return y, yc, x_[:, 1], s.simTime

    y = s.to_head(np.array(s.interpolateNN(x_)))
    yc = rich_cyls[-1].getSolutionHead()
    return yc, y, x_[:, 1], s.simTime


if __name__ == "__main__":

    N = 41

    sand = [0.045, 0.43, 0.15, 3, 1000, "Sand"]
    loam = [0.08, 0.43, 0.04, 1.6, 50, "Loam"]
    clay = [0.1, 0.4, 0.01, 1.1, 10, "Clay"]

    sim_times = np.linspace(0, 25, 250)  # temporal resolution of 0.1 d

    if rank == 0:
        fig, ax = plt.subplots(2, 3, figsize = (14, 14))
        t0 = timeit.default_timer()

    jobs = ([sand, 0.1, 0, 0], [loam, 0.1, 0, 1], [clay, 0.1, 0, 2], [sand, 0.05, 1, 0], [loam, 0.05, 1, 1], [clay, 0.05, 1, 2])

    for soil, qj, i, j in jobs:
        y, yc, x, t = solve(soil, sim_times, qj, N)
        if rank == 0:
            ax[i, j].plot(x, y, "r*", label = "sink based")
            ax[i, j].plot(np.linspace(r_root, r_out, N), yc, "b", label = "cylindrical")
            ax[i, j].set_xlabel("r (cm)")
            ax[i, j].set_ylabel("water potential (cm)")
            ax[i, j].title.set_text(soil[5] + ", q = {:g} cm/d, final: {:g} d".format(qj, t))
            ax[i, j].set_ylim(-15000, 0.)
            ax[i, j].legend()

    if rank == 0:
        print("Elapsed time: ", timeit.default_timer() - t0, "s")
        plt.show()

