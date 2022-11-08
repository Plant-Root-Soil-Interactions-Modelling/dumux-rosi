""" 

Functions to simplify setup of the scenarios for INARI


dumux NC story: 

    enum BCTypes {
        constantPressure = 1,
        constantConcentration = 1,
        constantFlux = 2,
        constantFluxCyl = 3,
        atmospheric = 4,
        freeDrainage = 5,
        outflow = 6,
        linear = 7,
        michaelisMenten = 8
    };



"""

import sys; sys.path.append("../../build-cmake/cpp/python_binding/"); sys.path.append("../modules/");
from hamcrest.core.core.isnone import none
sys.path.append("../../../CPlantBox/src/python_modules"); sys.path.append("../../../CPlantBox/");

import numpy as np
import timeit
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()

import plantbox as pb  # CPlantBox
import van_genuchten as vg
from xylem_flux import *

from rosi_richardsnc import RichardsNCSP  # C++ part (Dumux binding), macroscopic soil model
from rosi_richards import RichardsSP  # C++ part (Dumux binding), macroscopic soil model
from richards import RichardsWrapper  # Python part, macroscopic soil model


def create_soil_model(soil_, min_b , max_b , cell_number, p_top, p_bot, type, times = None, net_inf = None):
    """
        Creates a soil domain from @param min_b to @param max_b with resolution @param cell_number
        soil type is fixed and homogeneous 
        domain is periodic (if 2d or 3d)
        initial potentials are linear from @param p_top to @param p_bot
        
        returns soil_model (RichardsWrapper(RichardsSP())) and soil parameter (vg.Parameters)
    """
    soil = vg.Parameters(soil_)
    vg.create_mfp_lookup(soil, -1.e5, 1000)

    if type == 1:
        s = RichardsWrapper(RichardsSP())  # water only
    elif type == 2:
        s = RichardsWrapper(RichardsNCSP())  # water and one solute
    else:
        print("choose type, 1 = Richards, 2 = RichardsNCSP")

    s.initialize()
    s.createGrid(min_b, max_b, cell_number, True)  # [cm]

    # Initial conditions
    if type == 2:  # solute IC
        z_ = [0., -80., -80., -200.]
        v_ = [2.e-5, 2.e-5, 1.e-5, 1.e-5]  # [2.e-4, 2.e-4, 1.e-4, 1.e-4]  # TODO [0., 0., 0., 0.]  #
        s.setICZ_solute(v_[::-1], z_[::-1])  # ascending order...
    s.setLinearIC(p_top, p_bot)  # cm pressure head, equilibrium

    # BC
    if times:
        s.setTopBC("atmospheric", 0.5, [times, net_inf])  # 0.5 is dummy value
    else:
        s.setTopBC("noFlux")
    s.setBotBC("freeDrainage")

    if type == 2:  # solute BC
        sol_times = [0., 1., 1., 2., 2., 5., 5., 6., 6., 1.e3]
        sol_influx = -np.array([0., 0., 1.e-5, 1.e-5, 0., 0., 1.e-5, 1.e-5, 0., 0.])
        s.setTopBC_solute("managed", 0.5, [sol_times, sol_influx])
        # s.setTopBC_solute("constantFlux", 0.)
        s.setBotBC_solute("outflow", 0.)

    s.setVGParameters([soil_])
    s.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")

    # Solutes
    s.setParameter("Component.MolarMass", "6.2e-2")  # TODO no idea, where this is neeeded, i don't want to use moles ever (nitrate 62,0049 g/mol)
    s.setParameter("Component.LiquidDiffusionCoefficient", "1.7e-9")  # m2 s-1 # nitrate = 1700 um^2/sec
    # s.setParameter("Component.BufferPower", "140") # amonium has around 50?, no buffering for nitrate (in Kirk & Kronzuchiker 2005)

    s.initializeProblem()
    wilting_point = -15000
    s.setCriticalPressure(wilting_point)  # for boundary conditions constantFlow, constantFlowCyl, and atmospheric
    s.ddt = 1.e-5  # [day] initial Dumux time step

    return s, soil


def init_conductivities_const(r, kr_const = 1.8e-4, kx_const = 0.1):
    """ Hydraulic conductivities  kr [1/day], kx [cm3/day] """
    r.setKr([0, kr_const, kr_const, kr_const, kr_const, kr_const])
    r.setKx([1.e3, kx_const, kx_const, kx_const, kx_const, kx_const])


def init_conductivities_const_growth(r, kr_const = 1.8e-4, kx_const = 0.1):
    """ Hydraulic conductivities kr [1/day], kx [cm3/day] """
    kr = np.array([[-1e4, 0.], [-0.1, 0.], [0., kr_const], [1.e4, kr_const]])
    kx = np.array([[0, kx_const], [1e4, kx_const]])
    kr00 = np.array([[0., 0.]])  # artificial shoot
    kx00 = np.array([[0., 1.e3]])  # artificial shoot
    r.setKrTables([kr00[:, 1], kr[:, 1], kr[:, 1], kr[:, 1], kr[:, 1], kr[:, 1]],
                  [kr00[:, 0], kr[:, 0], kr[:, 0], kr[:, 0], kr[:, 0], kr[:, 0]])  # for each subtype
    r.setKxTables([kx00[:, 1], kx[:, 1], kx[:, 1], kx[:, 1], kx[:, 1], kx[:, 1]],
                  [kx00[:, 0], kx[:, 0], kx[:, 0], kx[:, 0], kx[:, 0], kx[:, 0]])  # for each subtype


def init_dynamic_simple_growth(r, kr0, kr1, kx0, kx1):
    """ a simplified parametrisation, based on init_dynamic_conductivities_growth """
    kr00 = np.array([[0., 0.]])  # artificial shoot
    kx00 = np.array([[0., 1.e3]])  # artificial shoot
    dt0 = 30.
    dt1 = 15.
    kr_f = 0.5
    kx_f = 5.
    kr0 = np.array([[-1e4, 0.], [-0.1, 0.], [0., kr0], [dt0, kr_f * kr0]])  # primals
    kr1 = np.array([[-1e4, 0.], [-0.1, 0.], [0., kr1], [dt1, kr_f * kr1]])  # laterals
    kx0 = np.array([[0., kx0], [dt0, kx_f * kx0]])  # primals
    kx1 = np.array([[0., kx1], [dt1, kx_f * kx1]])  # laterals
    r.setKrTables([kr00[:, 1], kr0[:, 1], kr1[:, 1], kr1[:, 1], kr0[:, 1], kr0[:, 1]],
                  [kr00[:, 0], kr0[:, 0], kr1[:, 0], kr1[:, 0], kr0[:, 0], kr0[:, 0]])
    r.setKxTables([kx00[:, 1], kx0[:, 1], kx1[:, 1], kx1[:, 1], kx0[:, 1], kx0[:, 1]],
                  [kx00[:, 0], kx0[:, 0], kx1[:, 0], kx1[:, 0], kx0[:, 0], kx0[:, 0]])


def init_maize_conductivities(r):
    """ Hydraulic conductivities for maize following Couvreur et al. (2012) originally from Doussan et al. (1998) """
    kr00 = np.array([[0., 0.]])  # artificial shoot
    kx00 = np.array([[0., 1.e3]])  # artificial shoot
    kr0 = np.array([[-1e4, 0.], [-0.1, 0.], [0., 0.000181], [8., 0.000181], [10, 0.0000648], [18, 0.0000648], [25, 0.0000173], [300, 0.0000173]])
    kr1 = np.array([[-1e4, 0.], [-0.1, 0.], [0., 0.000181], [10., 0.000181], [16, 0.0000173], [300, 0.0000173]])
    kx0 = np.array([[0., 0.0000864], [5., 0.00173], [12., 0.0295], [15., 0.0295], [20., 0.432], [300., 0.432]])
    kx1 = np.array([[0., 0.0000864], [5., 0.0000864], [10., 0.0000864], [12., 0.0006048], [20., 0.0006048], [23., 0.00173], [300., 0.00173]])
    r.setKrTables([kr00[:, 1], kr0[:, 1], kr1[:, 1], kr1[:, 1], kr0[:, 1], kr0[:, 1]],
                  [kr00[:, 0], kr0[:, 0], kr1[:, 0], kr1[:, 0], kr0[:, 0], kr0[:, 0]])
    r.setKxTables([kx00[:, 1], kx0[:, 1], kx1[:, 1], kx1[:, 1], kx0[:, 1], kx0[:, 1]],
                  [kx00[:, 0], kx0[:, 0], kx1[:, 0], kx1[:, 0], kx0[:, 0], kx0[:, 0]])


def init_lupine_conductivities(r):
    """ Hydraulic conductivities for lupine following Zarebanadkouki et al. (2016) """
    kr00 = np.array([[0., 0.]])  # artificial shoot
    kx00 = np.array([[0., 1.e3]])  # artificial shoot
    kr0 = np.array([[-1.e4, 0.], [-0.1, 0.], [0., 1.14e-03], [2, 1.09e-03], [4, 1.03e-03], [6, 9.83e-04], [8, 9.35e-04], [10, 8.90e-04],
                    [12, 8.47e-04], [14, 8.06e-04], [16, 7.67e-04], [18, 7.30e-04], [20, 6.95e-04], [22, 6.62e-04], [24, 6.30e-04], [26, 5.99e-04],
                    [28, 5.70e-04], [30, 5.43e-04], [32, 5.17e-04]])
    kr1 = np.array([[-1e4, 0.], [-0.1, 0.], [0., 4.11e-03], [1, 3.89e-03], [2, 3.67e-03], [3, 3.47e-03], [4, 3.28e-03],
                    [5, 3.10e-03], [6, 2.93e-03], [7, 2.77e-03], [8, 2.62e-03], [9, 2.48e-03], [10, 2.34e-03], [11, 2.21e-03],
                    [12, 2.09e-03], [13, 1.98e-03], [14, 1.87e-03], [15, 1.77e-03], [16, 1.67e-03], [17, 1.58e-03]])
    kx0 = np.array([[0., 6.74e-02], [2, 7.48e-02], [4, 8.30e-02], [6, 9.21e-02], [8, 1.02e-01], [10, 1.13e-01],
                    [12, 1.26e-01], [14, 1.40e-01], [16, 1.55e-01], [18, 1.72e-01], [20, 1.91e-01], [22, 2.12e-01], [24, 2.35e-01],
                    [26, 2.61e-01], [28, 2.90e-01], [30, 3.21e-01], [32, 3.57e-01]])
    kx1 = np.array([[0., 4.07e-04], [1, 5.00e-04], [2, 6.15e-04], [3, 7.56e-04], [4, 9.30e-04], [5, 1.14e-03],
                    [6, 1.41e-03], [7, 1.73e-03], [8, 2.12e-03], [9, 2.61e-03], [10, 3.21e-03], [11, 3.95e-03], [12, 4.86e-03],
                    [13, 5.97e-03], [14, 7.34e-03], [15, 9.03e-03], [16, 1.11e-02], [17, 1.36e-02]])
    r.setKrTables([kr00[:, 1], kr0[:, 1], kr1[:, 1], kr1[:, 1], kr0[:, 1], kr0[:, 1]],
                  [kr00[:, 0], kr0[:, 0], kr1[:, 0], kr1[:, 0], kr0[:, 0], kr0[:, 0]])
    r.setKxTables([kx00[:, 1], kx0[:, 1], kx1[:, 1], kx1[:, 1], kx0[:, 1], kx0[:, 1]],
                  [kx00[:, 0], kx0[:, 0], kx1[:, 0], kx1[:, 0], kx0[:, 0], kx0[:, 0]])


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
    return XylemFluxPython(rs)


def set_all_sd(rs, s):
    """ # sets all standard deviation to a percantage, i.e. value*s """
    for p in rs.getRootRandomParameter():
        p.a_s = p.a * s
        p.lbs = p.lb * s
        p.las = p.la * s
        p.lns = p.ln * s
        p.lmaxs = p.lmax * s
        p.rs = p.r * s
        p.thetas = p.theta * s
        p.rlts = p.rlt * s  # no used
        p.ldelays = p.ldelay * s
    seed = rs.getRootSystemParameter()  # SeedRandomParameter
    seed.firstBs = seed.firstB * s
    seed.delayBs = seed.delayB * s
    seed.maxBs = seed.maxB * s
    seed.firstSBs = seed.firstSB * s
    seed.delaySBs = seed.delaySB * s
    seed.delayRCs = seed.delayRC * s
    seed.nCs = seed.nCs * s
    seed.nzs = seed.nzs * s
    # todo seed position s


def create_mapped_rootsystem(min_b , max_b , cell_number, soil_model, fname, stochastic = False):
    """ loads a rmsl file, or creates a rootsystem opening an xml parameter set,  
        and maps it to the soil_model """
    global picker  # make sure it is not garbage collected away...

    if fname.endswith(".rsml"):
        r = XylemFluxPython(fname)
    elif fname.endswith(".xml"):
        if rank == 0:
            if stochastic:
                seed = np.random.randint(0, 1e6)
            else:
                seed = 1  # always the same random seed
        else:
            seed = None
        seed = comm.bcast(seed, root = 0)  # random seed must be the same for each process
        rs = pb.MappedRootSystem()
        rs.setSeed(seed)
        rs.readParameters(fname)
        if not stochastic:
            set_all_sd(rs, 0.)
        rs.initializeDB(4, 5)
        rs.simulate(1., True)
        rs.setGeometry(pb.SDF_PlantBox(1.e6, 1.e6, np.abs(min_b[2])))
        r = XylemFluxPython(rs)

        # print("HERE***********************************")
        # print([s.x for s in r.rs.segments])
        # print([s.y for s in r.rs.segments])
        # # for i in range(0, len(r.rs.segments)): # ????????
        # #     print(r.rs.seg2cell[i])

    r.rs.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                            pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), cut = False)

    # print("HERE***********************************")
    # print([s.x for s in r.rs.segments])
    # print([s.y for s in r.rs.segments])
    # ss
    comm.barrier()
    # print("survived setRectangularGrid", rank)

    picker = lambda x, y, z: soil_model.pick([x, y, z])  #  function that return the index of a given position in the soil grid (should work for any grid - needs testing)
    r.rs.setSoilGrid(picker)  # maps segments, maps root segements and soil grid indices to each other in both directions
    comm.barrier()
    # print("survived setSoilGrid", rank)

    if rank == 0:
        init_conductivities_const(r)
    return r


def write_files(file_name, psi_x, psi_i, sink, times, trans, psi_s, conc = None):
    """  saves numpy arrays ass npy files """
    np.save('results/psix_' + file_name, np.array(psi_x))  # xylem pressure head per segment [cm]
    np.save('results/psiinterface_' + file_name, np.array(psi_i))  # pressure head at interface per segment [cm]
    np.save('results/sink_' + file_name, -np.array(sink))  # sink per segment [cm3/day]
    np.save('results/transpiration_' + file_name, np.vstack((times, -np.array(trans))))  # time [day], transpiration [cm3/day]
    np.save('results/soil_' + file_name, np.array(psi_s))  # soil potential per cell [cm]
    if conc is not None:
        np.save('results/soilc_' + file_name, np.array(conc))  # soil potential per cell [cm]


def simulate_const(s, r, trans, sim_time, dt):
    """ 
        classic model:
        potential at root soil interface equals mean matric potential of surrounding finite volume
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
