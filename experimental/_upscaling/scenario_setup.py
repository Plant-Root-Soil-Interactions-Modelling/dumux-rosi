""" functions to simplify setup of the scenarios """

import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

import numpy as np
import timeit
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()

import plantbox as pb  # CPlantBox
import functional.van_genuchten as vg
from functional.xylem_flux import sinusoidal2
from ups import XylemFluxUps  # == XylemFluxPython with upscaling stuff added

from rosi_richards import RichardsSP  # C++ part (Dumux binding), macroscopic soil model
from richards import RichardsWrapper  # Python part, macroscopic soil model


def create_soil_model(soil_, min_b , max_b , cell_number, p_top, p_bot):
    """
        Creates a soil domain from @param min_b to @param max_b with resolution @param cell_number
        soil type is fixed and homogeneous 
        domain is periodic (if 2d or 3d)
        initial potentials are linear from @param p_top to @param p_bot
        
        returns soil_model (RichardsWrapper(RichardsSP())) and soil parameter (vg.Parameters)
    """
    soil = vg.Parameters(soil_)
    vg.create_mfp_lookup(soil, -1.e5, 1000)
    # if (cell_number[0] == 1 and cell_number[1] == 1):  # 1D
    #     periodic = False
    # else:  # 2D or 3D
    #     periodic = True
    periodic = True
    s = RichardsWrapper(RichardsSP())
    s.initialize()
    s.createGrid(min_b, max_b, cell_number, periodic)  # [cm]
    s.setLinearIC(p_top, p_bot)  # cm pressure head, equilibrium
    s.setTopBC("noFlux")
    s.setBotBC("freeDrainage")
    s.setVGParameters([soil_])
    s.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
    s.initializeProblem()
    wilting_point = -15000
    s.setCriticalPressure(wilting_point)  # for boundary conditions constantFlow, constantFlowCyl, and atmospheric
    s.ddt = 1.e-5  # [day] initial Dumux time step

    return s, soil


def init_conductivities_const(r):
    """ Hydraulic conductivities - for Jans scenarios, but constant """
    # Scenario 1
    kr_const = 1.8e-4  # [1/day]
    kx_const = 0.1  # [cm3/day]
    r.setKr([kr_const])
    r.setKx([kx_const])


def create_mapped_singleroot(min_b , max_b , cell_number, soil_model, ns = 100, l = 50 , a = 0.05):
    """ creates a single root mapped to a soil with @param ns segments, length l, and radius a """
    global picker  # make sure it is not garbage collected away...
    r = create_singleroot(ns, l, a)
    r.rs.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                            pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), cut = False)
    picker = lambda x, y, z: soil_model.pick([x, y, z])  #  function that return the index of a given position in the soil grid (should work for any grid - needs testing)
    r.rs.setSoilGrid(picker)  # maps segments, maps root segements and soil grid indices to each other in both directions
    init_conductivities_const(r)
    return r


def create_singleroot(ns = 100, l = 50 , a = 0.05):
    """ creates a single root with @param ns segments, length l, and radius a """
    radii = np.array([a] * ns)
    nodes = [pb.Vector3d(0, 0, 0)]
    segs = []
    dx = l / ns
    z_ = np.linspace(-dx, -l , ns)
    for i in range(0, ns):
        nodes.append(pb.Vector3d(0, 0, z_[i]))
        segs.append(pb.Vector2i(i, i + 1))
    rs = pb.MappedSegments(nodes, segs, radii)
    return XylemFluxUps(rs)


def create_mapped_rootsystem(min_b , max_b , cell_number, soil_model, fname):
    """ loads a rsml file and maps it to the soil_model """
    global picker  # make sure it is not garbage collected away...
    r = XylemFluxUps(fname)  # see rootsystem.py (in upscaling)
    r.rs.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                            pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), cut = False)
    picker = lambda x, y, z: soil_model.pick([x, y, z])  #  function that return the index of a given position in the soil grid (should work for any grid - needs testing)
    r.rs.setSoilGrid(picker)  # maps segments, maps root segements and soil grid indices to each other in both directions
    init_conductivities_const(r)
    return r


def write_files(file_name, psi_x, psi_i, sink, times, trans, psi_s):
    """  saves numpy arrays ass npy files """
    np.save('results/psix_' + file_name, np.array(psi_x))  # xylem pressure head per segment [cm]
    np.save('results/psiinterface_' + file_name, np.array(psi_i))  # pressure head at interface per segment [cm]
    np.save('results/sink_' + file_name, -np.array(sink))  # sink per segment [cm3/day]
    np.save('results/transpiration_' + file_name, np.vstack((times, -np.array(trans))))  # time [day], transpiration [cm3/day]
    np.save('results/soil_' + file_name, np.array(psi_s))  # soil potential per cell [cm]


def simulate_const(s, r, trans, sim_time, dt):
    """ 
    TODO 
    """
    wilting_point = -15000  # cm
    skip = 6  # for output and results, skip iteration
    rs_age = 0.  # day

    start_time = timeit.default_timer()
    psi_x_, psi_s_, sink_ , x_, y_, psi_s2_ = [], [], [], [], [], []  # for post processing
    sx = s.getSolutionHead()  # inital condition, solverbase.py
    ns = len(r.rs.segments)
    if rank == 0:
        map_ = r.rs.seg2cell  # because seg2cell is a map
        mapping = np.array([map_[j] for j in range(0, ns)], dtype = np.int64)  # convert to a list

    N = int(np.ceil(sim_time / dt))

    """ simulation loop """
    for i in range(0, N):

        t = i * dt  # current simulation time

        """ 1. xylem model """
        if rank == 0:  # Root part is not parallel
            rx = r.solve(rs_age, -trans * sinusoidal2(t, dt), 0., sx, cells = True, wilting_point = wilting_point)  # xylem_flux.py
            fluxes = r.soilFluxes(rs_age, rx, sx, False)  # class XylemFlux is defined in MappedOrganism.h, approx = False
        else:
            fluxes = None
        fluxes = comm.bcast(fluxes, root = 0)  # Soil part runs parallel

        """ 2. soil model """
        s.setSource(fluxes.copy())  # richards.py
        s.solve(dt)
        sx = s.getSolutionHead()  # richards.py

        """ validity check """

        """ remember results ... """
        if rank == 0 and i % skip == 0:

            sx_ = sx[:, 0]
            psi_x_.append(rx.copy())  # cm (per root node)
            psi_s_.append(np.array([sx_[ci] for ci in mapping]))  # cm (per root segment)
            sink = np.zeros(sx_.shape)
            for k, v in fluxes.items():
                sink[k] += v
            sink_.append(sink)  # cm3/day (per soil cell)
            x_.append(t)  # day
            y_.append(np.sum(sink))  # cm3/day
            psi_s2_.append(sx_)  # cm (per soil cell)

            min_sx, min_rx, max_sx, max_rx = np.min(sx), np.min(rx), np.max(sx), np.max(rx)
            n = round(float(i) / float(N) * 100.)
            print("\n[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], [{:g}, {:g}] cm soil [{:g}, {:g}] cm root at {:g}, {:g}"
                    .format(min_sx, max_sx, min_rx, max_rx, np.sum(sink), -trans * sinusoidal2(t, dt)))

    if rank == 0:
        print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")

    return psi_x_, psi_s_, sink_, x_, y_, psi_s2_


if __name__ == '__main__':

    theta_r = 0.025
    theta_s = 0.403
    alpha = 0.0383  # (cm-1) soil
    n = 1.3774
    k_sat = 60.  # (cm d-1)
    soil_ = [theta_r, theta_s, alpha, n, k_sat]
    s, soil = create_soil_model(soil_, [-1, -1, -150.], [1, 1, 0.], [1, 1, 55], -310, -200)

    print()
    print(s)
    print(soil)

    """ TODO: tests would be nice, or a minimal example setup ... """
