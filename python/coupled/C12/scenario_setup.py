""" 
Basic set up for Benchmark C12 for a static root system
"""
import plantbox as pb
from functional.root_conductivities import init_conductivities

import functional.van_genuchten as vg
from functional.xylem_flux import XylemFluxPython  # Python hybrid solver
from functional.HydraulicsDoussan import HydraulicsDoussan  # Doussan solver

from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part

import numpy as np
from scipy.interpolate import RegularGridInterpolator


def get_domain1D():
    min_b = [-4., -4., -15.]
    max_b = [4., 4., 0.]
    cell_number = [1, 1, 15]
    return min_b, max_b, cell_number


def get_domain3D():
    min_b = [-4., -4., -15.]
    max_b = [4., 4., 0.]
    cell_number = [7, 7, 15]
    return min_b, max_b, cell_number


def set_scenario(dimension, age_dependent = False):

    if dimension == "1D":
        min_b, max_b, cell_number = get_domain1D()
        slope = "100"  # big leads to strange behaviour
    elif dimension == "3D":
        min_b, max_b, cell_number = get_domain3D()
        slope = "500"
    else:
        raise

    initial = -659.8 - min_b[2] / 2  # + 7.5

    loam = [0.08, 0.43, 0.04, 1.6, 50]
    soil = loam
    sra_table_lookup = open_sra_lookup("../table_loam")

    wilting_point = -15000  # cm

    """ Initialize macroscopic soil model """
    sp = vg.Parameters(soil)  # for debugging
    s = RichardsWrapper(RichardsSP())
    s.initialize()
    s.createGrid(min_b, max_b, cell_number, False)  # [cm]
    s.setHomogeneousIC(initial, True)  # cm pressure head, equilibrium
    s.setTopBC("noFlux")
    s.setBotBC("noFlux")
    s.setVGParameters([soil])
    s.setParameter("Newton.EnableChop", "True")
    s.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
    s.setParameter("Soil.SourceSlope", slope)  # turns regularisation of the source term on, will change the shape of actual transpiration...
    s.initializeProblem()
    s.setCriticalPressure(wilting_point)
    s.ddt = 1.e-5  # [day] initial Dumux time step

    """ root hydraulic model"""
    fname = "../../../grids/RootSystem8.rsml"
    trans = 6.4  # cm3 /day (sinusoidal)
    r = HydraulicsDoussan(fname)
    r.rs.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                            pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), True)  # cutting
    init_conductivities(r, age_dependent)
    rs_age = np.max(r.get_ages())

    """ coupling (map indices) """
    if dimension == "1D":
       picker = lambda x, y, z: s.pick([0., 0., z])
    elif dimension == "3D":
        picker = lambda x, y, z: s.pick([x, y, z])

    r.rs.setSoilGrid(picker)  # maps segment
    seg2cell = r.rs.seg2cell
    ns = len(r.rs.segments)
    mapping = np.array([seg2cell[j] for j in range(0, ns)])

    """ simulation time """
    sim_time = 7.1  # 0.65  # 0.25  # [day]
    dt = 360 / (24 * 3600)  # time step [day]
    skip = 1  # for output and results, skip iteration

    return r, rs_age, trans, wilting_point, soil, s, sra_table_lookup, mapping, sim_time, dt, skip


def sinusoidal(t):
    return np.sin(2. * np.pi * np.array(t) - 0.5 * np.pi) + 1.


def soil_root_interface_table(rx, sx, inner_kr_, rho_, f):
    """
    finds potential at the soil root interface
        
    rx             xylem matric potential [cm]
    sx             bulk soil matric potential [cm]
    inner_kr       root radius times hydraulic conductivity [cm/day] 
    rho            geometry factor [1]
    f              function to look up the potentials
    """
    try:
        rsx = f((rx, sx, inner_kr_ , rho_))
    except Exception as err:
        print("rx", rx.shape, np.min(rx), np.max(rx))  # 0, -16000
        print("sx", sx.shape, np.min(sx), np.max(sx))  # 0, -16000
        print("inner_kr_", inner_kr_.shape, np.min(inner_kr_), np.max(inner_kr_))  # 1.e-7 - 1.e-4
        print("rho", rho_.shape, np.min(rho_), np.max(rho_))  # 1. - 200.
        raise err

    return rsx


def open_sra_lookup(filename):
    """ opens the look from a file """
    sra_table = np.load(filename + ".npy")
    x = np.load(filename + "_.npy", allow_pickle = True)
    kx_ = x[0]
    sx_ = x[1]
    inner_ = x[2]
    outer_ = x[3]
    return RegularGridInterpolator((kx_, sx_, inner_, outer_), sra_table)

