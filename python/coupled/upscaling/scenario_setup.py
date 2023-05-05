""" 
Basic set-ups for the new upscaling manuscript (DL 29.3.2023)
"""
import plantbox as pb
from functional.root_conductivities import init_conductivities

import functional.van_genuchten as vg
from functional.xylem_flux import XylemFluxPython  # Python hybrid solver
from functional.HydraulicsDoussan import HydraulicsDoussan  # Doussan solver
from functional.Perirhizal import PerirhizalPython

from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part

import numpy as np
from scipy.interpolate import RegularGridInterpolator


def soil_vg_(name:str):
    """ 
    Van Genuchten parameter for soil from Hydrus1D, 
    called by maize() and soybean() 
    
    4D look up tables are created with teh script create_sra_table_v2
    """
    soil = {}
    soil["jan_comp"] = [0.025, 0.403, 0.0383, 1.3774, 60.]
    soil["hydrus_loam"] = [0.078, 0.43, 0.036, 1.56, 24.96]
    soil["hydrus_clay"] = [0.068, 0.38, 0.008, 1.09, 4.8]
    soil["hydrus_sand"] = [0.045, 0.43, 0.145, 2.68, 712.8]
    soil["hydrus_sandyloam"] = [0.065, 0.41, 0.075, 1.89, 106.1]
    table_name = "table_{:s}".format(name)  # name for 4D look up table ########################
    return soil[name], table_name


def maize_(dim:str):
    """ parameters for maize simulation """
    min_b = np.array([-38., -8., -100.])
    max_b = np.array([38., 8., 0.])
    if dim == "1D":
        cell_number = np.array([1, 1, 100])
    else:
        cell_number = np.array([76, 16, 100])
    return min_b, max_b, cell_number


def soybean_(dim:str):
    """ parameters for soybean simulation """
    min_b = np.array([-38, -1.5, -100.])
    max_b = np.array([38, 1.5, 0.])
    if dim == "1D":
        cell_number = np.array([1, 1, 100])
    else:
        cell_number = np.array([76, 3, 100])
    return min_b, max_b, cell_number


def maize_conductivities(r, skr = 1., skx = 1.):
    """ Hydraulic conductivities for maize modified from Couvreur et al. (2012), originally from Doussan et al. (1998) """

    kr00 = np.array([[0., 0.]])  # artificial shoot
    kx00 = np.array([[0., 1.e3]])  # artificial shoot

    kr0 = np.array([[-1e4, 0.], [-0.1, 0.], [0., 0.000181], [8., 0.000181], [10, 0.0000648], [18, 0.0000648], [25, 0.0000173], [300, 0.0000173]])
    kr1 = np.array([[-1e4, 0.], [-0.1, 0.], [0., 0.000181], [10., 0.000181], [16, 0.0000173], [300, 0.0000173]])

    kx0 = np.array([[0., 0.000864], [5., 0.00173], [12., 0.0295], [15., 0.0295], [20., 0.432], [300., 0.432]])
    kx1 = np.array([[0., 0.000864], [5., 0.0000864], [10., 0.0000864], [12., 0.0006048], [20., 0.0006048], [23., 0.00173], [300., 0.00173]])

    kr01 = np.minimum(skr * kr0[:, 1], 1.)  # kr0[:, 1] are values
    kr11 = np.minimum(skr * kr1[:, 1], 1.)  # kr1[:, 1] are values
    r.setKrTables([kr00[:, 1], kr01, kr11, kr11, kr01, kr01],
                  [kr00[:, 0], kr0[:, 0], kr1[:, 0], kr1[:, 0], kr0[:, 0], kr0[:, 0]])
    kx01 = np.minimum(skx * kx0[:, 1], 1.)
    kx11 = np.minimum(skx * kx1[:, 1], 1.)
    r.setKxTables([kx00[:, 1], kx01, kx11, kx11, kx01, kx01],
                  [kx00[:, 0], kx0[:, 0], kx1[:, 0], kx1[:, 0], kx0[:, 0], kx0[:, 0]])


def lupine_conductivities(r, skr = 1., skx = 1.):
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

    kr01 = np.minimum(skr * kr0[:, 1], 1.)
    kr11 = np.minimum(skr * kr1[:, 1], 1.)
    r.setKrTables([kr00[:, 1], kr01, kr11, kr11, kr01, kr01],
                  [kr00[:, 0], kr0[:, 0], kr1[:, 0], kr1[:, 0], kr0[:, 0], kr0[:, 0]])
    kx01 = np.minimum(skx * kx0[:, 1], 1.)
    kx11 = np.minimum(skx * kx1[:, 1], 1.)
    r.setKxTables([kx00[:, 1], kx01, kx11, kx11, kx01, kx01],
                  [kx00[:, 0], kx0[:, 0], kx1[:, 0], kx1[:, 0], kx0[:, 0], kx0[:, 0]])


def set_scenario(plant, dim, initial, soil, outer_method):
    """ 
    Sets up a Scenario     
    
    plant            plant name: 'maize' or 'soybean'
    dim              soil model dimensionality: '1D' or '3D'
    initial          initial soil total potential [cm] 
    soil             name of the soil 4D look up table, see soil_vg_() for the names and corresponding VG parameters
    outer_method     method to determine outer perirhizal radii ('voronoi', 'length', 'surface', or 'volume')    
    """
    assert plant == "maize" or plant == "soybean", "plant should be 'maize', or 'soybean' "
    assert dim == "3D" or dim == "1D", "dim should be '1D' or '3D'"
    assert soil in ["hydrus_loam", "hydrus_clay", "hydrus_sand", "hydrus_sandyloam"], "soil should be 'hydrus_loam', 'hydrus_clay', 'hydrus_sand' or 'hydrus_sandyloam' "
    assert outer_method in ["voronoi", "length", "surface", "volume"], "outer_method should be 'voronoi', 'length', 'surface', or 'volume'"

    wilting_point = -15000  # cm
    random_seed = 1  # random seed
    slope = "1000"  # cm
    trans_ = 0.5  # cm / day
    rs_age = 21.  # initial age in days

    soil_, table_name = soil_vg_(soil)
    if plant == "maize":
        min_b, max_b, cell_number = maize_(dim)
        param_name = "Zeamays_synMRI_modified.xml"
    elif plant == "soybean":
        min_b, max_b, cell_number = soybean_(dim)
        param_name = "Glycine_max_Moraes2020_opt2_modified.xml"

    trans = (max_b[0] - min_b[0]) * (max_b[1] - min_b[1]) * trans_  # cm3 / day
    initial -= min_b[2] / 2

    """ initialize macroscopic soil model """
    soil = vg.Parameters(soil_)
    vg.create_mfp_lookup(soil, -1.e5, 1000)
    sra_table_lookup = open_sra_lookup("../" + table_name)
    s = RichardsWrapper(RichardsSP())
    s.initialize()
    if dim == "1D":
        s.createGrid(min_b, max_b, cell_number, periodic = False)
    elif dim == "3D":
        s.createGrid(min_b, max_b, cell_number, periodic = True)
    s.setHomogeneousIC(initial, True)  # cm pressure head top, equilibrium (contant total potential)
    s.setTopBC("noFlux")
    s.setBotBC("noFlux")
    s.setVGParameters([soil_])
    s.setParameter("Newton.EnableChop", "True")
    s.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
    s.setParameter("Soil.SourceSlope", slope)
    s.initializeProblem()
    s.setCriticalPressure(wilting_point)
    s.ddt = 1.e-5  # [day] initial Dumux time step

    """ root hydraulic model"""
    rs = pb.MappedRootSystem()
    rs.setSeed(random_seed)
    rs.readParameters("data/" + param_name)
    rs.setGeometry(pb.SDF_PlantBox(1.e6, 1.e6, np.abs(min_b[2])))
    # rs.initializeDB(4, 5)
    rs.initialize()
    rs.simulate(rs_age, True)

    r = HydraulicsDoussan(rs)

    r.rs.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                            pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), False)  # cutting
    if plant == "maize":
        maize_conductivities(r, 1., 1.)
    elif plant == "soybean":
        lupine_conductivities(r, 1., 1.)

    """ coupling roots to macroscopic soil """
    if dim == "1D":
        picker = lambda x, y, z: s.pick([0., 0., z])
    elif dim == "3D":
        picker = lambda x, y, z: s.pick([x, y, z])

    r.rs.setSoilGrid(picker)  # maps segment
    seg2cell = r.rs.seg2cell
    ns = len(r.rs.segments)
    mapping = np.array([seg2cell[j] for j in range(0, ns)])

    """ outer radii """
    if outer_method == "voronoi":
        outer_ = PerirhizalPython(rs).get_outer_radii_voronoi()
        print("Voroni open or outside regions", np.count_nonzero(np.isnan(outer_)), "/", len(outer_), "are replaced by mean value")
        outer_mean = np.nanmean(outer_) * np.ones(outer_.shape)
        outer_[np.isnan(outer_)] = outer_mean[np.isnan(outer_)]
        outer_ = outer_[1:]  # nodes to segs
    else:
        outer_ = PerirhizalPython(rs).get_outer_radii(outer_method)

    inner_ = rs.radii
    rho = np.divide(outer_, np.array(inner_))
    rho = np.expand_dims(rho, axis = 1)

    return r, rho, rs_age, trans, wilting_point, soil, s, sra_table_lookup, mapping


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
        print("An exception occured in soil_root_interface_table() with values:")
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


def write_files(file_name, hx, hsr, sink, times, trans, trans2, hs):
    """  saves numpy arrays as npy files """
    np.save('results/hx_' + file_name, np.array(hx))  # xylem pressure head per segment [cm]
    np.save('results/hsr_' + file_name, np.array(hsr))  # pressure head at interface per segment [cm]
    np.save('results/sink_' + file_name, -np.array(sink))  # sink per soil cell [cm3/day]
    np.save('results/transpiration_' + file_name, np.vstack((times, -np.array(trans), -np.array(trans2))))  # time [day], transpiration [cm3/day]
    np.save('results/hs_' + file_name, np.array(hs))  # soil water matric potential per soil cell [cm]

