import sys; sys.path.append("../../python/modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

import functional.van_genuchten as vg
from functional.xylem_flux import XylemFluxPython  # Python hybrid solver
from functional.HydraulicsDoussan import HydraulicsDoussan  # Doussan solver
from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part
import plantbox as pb
import aggregated_rs as agg

import numpy as np
from scipy.interpolate import RegularGridInterpolator

"""




"""


def get_domain1D():
    min_b = [-7.5, -37.5 / 2, -110.]
    max_b = [7.5, 37.5 / 2, 0.]
    cell_number = [1, 1, 55]  # [8, 38, 55]  # 2cm3
    return min_b, max_b, cell_number


def set_scenario1D(sstr, doussan = False):

    """ soil """
    min_b, max_b, cell_number = get_domain1D()

    if sstr == "_wet":
        p_top = -310
        p_bot = -200
    elif sstr == "_dry":
        p_top = -5000
        p_bot = -200
    else:
        raise "Unknown scenario souldbe '_wet' or '_dry'"

    alpha = 0.018;  # (cm-1)
    n = 1.8;
    Ks = 28.46;  # (cm d-1)
    loam = [0.08, 0.43, alpha, n, Ks]
    soil_ = loam
    soil = vg.Parameters(soil_)
    sra_table_lookup = open_sra_lookup("../../python/coupled/table_jan2")

    """ Initialize macroscopic soil model """
    s = RichardsWrapper(RichardsSP())
    s.initialize()
    s.createGrid(min_b, max_b, cell_number, False)  # [cm]
    s.setLinearIC(p_top, p_bot)  # cm pressure head, equilibrium
    s.setTopBC("noFlux")
    s.setBotBC("noFlux")
    s.setVGParameters([soil_])
    s.setParameter("Newton.EnableChop", "True")
    s.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
    # s.setParameter("Soil.SourceSlope", "1000")
    s.initializeProblem()
    wilting_point = -15000  # [cm]
    s.setCriticalPressure(wilting_point)
    s.ddt = 1.e-5  # [day] initial Dumux time step

    """ root system """
    fname = "../../grids/RootSystem_verysimple2.rsml"
    trans = 0.5 * 15 * 75  # average per day [cm3 /day] (sinusoidal) ########################################################### ANPASSEN
    rs_age = 78  # initial root system age

    """ root hydraulic model"""
    if doussan:
        r = HydraulicsDoussan(fname)
    else:
        r = XylemFluxPython(fname)

    types = r.rs.subTypes  # simplify root types
    types = (np.array(types) >= 12) * 1  # all roots type 0, only >=12 are laterals type 1
    r.rs.subTypes = list(types)
    r.rs.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                            pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), True)  # cutting
    agg.init_conductivities_const(r)

    """ coupling (map indices) """
    picker = lambda x, y, z: s.pick([0., 0., z])
    r.rs.setSoilGrid(picker)  # maps segment
    seg2cell = r.rs.seg2cell
    ns = len(r.rs.segments)
    mapping = np.array([seg2cell[j] for j in range(0, ns)])

    """ simulation time """
    sim_time = 7.1  # 0.65  # 0.25  # [day]
    dt = 360 / (24 * 3600)  # time step [day]
    skip = 1  # for output and results, skip iteration

    return r, rs_age, trans, wilting_point, soil_, s, sra_table_lookup, mapping, sim_time, dt, skip


def sinusoidal(t):
    return np.sin(2. * np.pi * np.array(t) - 0.5 * np.pi) + 1.


def soil_root_interface_table2(rx, sx, inner_kr_, rho_, f):
    assert rx.shape == sx.shape
    rsx = f((rx, sx, inner_kr_ , rho_))
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

