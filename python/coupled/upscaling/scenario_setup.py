""" 
Basic set-ups for the new upscaling manuscript (DL 29.3.2023)
"""
import plantbox as pb
from functional.root_conductivities import init_conductivities
import functional.van_genuchten as vg
from functional.xylem_flux import XylemFluxPython  # Python hybrid solver
from functional.Perirhizal import PerirhizalPython

from functional.PlantHydraulicParameters import PlantHydraulicParameters  # Doussan solver
from functional.PlantHydraulicModel import PlantHydraulicModel  # Doussan solver

from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part

from conductivities import *

import sys
import numpy as np
from scipy.interpolate import RegularGridInterpolator


def sinusoidal(t):
    """ sinusoidal function (used for transpiration) (integral over one day is 1)"""
    return np.sin(2. * np.pi * np.array(t) - 0.5 * np.pi) + 1.


def sinusoidal2(t, dt):
    """ sinusoidal function from 6:00 - 18:00, 0 otherwise (integral over one day is 1)"""
    return np.maximum(0., np.pi * (np.cos(2 * np.pi * (t - 0.5)) + np.cos(2 * np.pi * ((t + dt) - 0.5))) / 2)


def soil_vg_(name:str):
    """ 
    Van Genuchten parameter for soil from Hydrus1D, 
    called by maize() and soybean() 
    
    4D look up tables are created with thes script create_sra_table_v2
    """
    soil = {}
    soil["jan_comp"] = [0.025, 0.403, 0.0383, 1.3774, 60.]
    soil["hydrus_loam"] = [0.078, 0.43, 0.036, 1.56, 24.96]
    soil["hydrus_clay"] = [0.068, 0.38, 0.008, 1.09, 4.8]
    soil["hydrus_sand"] = [0.045, 0.43, 0.145, 2.68, 712.8]
    soil["hydrus_sandyloam"] = [0.065, 0.41, 0.075, 1.89, 106.1]
    table_name = "table_{:s}".format(name)  # name for 4D look up table ########################
    return soil[name], table_name


def springbarley_(dim:str):
    """ parameters for maize simulation """
    min_b = np.array([-6.5, -1.5, -150.])
    max_b = np.array([6.5, 1.5, 0.])
    if dim == "1D":
        cell_number = np.array([1, 1, 150])
    elif dim == "2D":
        cell_number = np.array([13, 1, 150])
    else:
        cell_number = np.array([13, 3, 150])
    return min_b, max_b, cell_number


def maize_(dim:str):
    """ parameters for maize simulation """
    min_b = np.array([-38., -8., -150.])
    max_b = np.array([38., 8., 0.])
    if dim == "1D":
        cell_number = np.array([1, 1, 150])
    elif dim == "2D":
        cell_number = np.array([38, 1, 150])
    else:
        cell_number = np.array([76, 16, 150])
    return min_b, max_b, cell_number


def soybean_(dim:str):
    """ parameters for soybean simulation """
    min_b = np.array([-38, -1.5, -150.])
    max_b = np.array([38, 1.5, 0.])
    if dim == "1D":
        cell_number = np.array([1, 1, 150])
    else:
        cell_number = np.array([76, 3, 150])
    return min_b, max_b, cell_number


def set_scenario(plant, dim, initial, soil, outer_method):
    """ 
    Sets up a Scenario     
    
    plant            plant name: 'maize' or 'soybean' or 'springbarley'
    dim              soil model dimensionality: '1D' or '3D'
    initial          initial soil total potential [cm] 
    soil             name of the soil 4D look up table, see soil_vg_() for the names and corresponding VG parameters
    outer_method     method to determine outer perirhizal radii ('voronoi', 'length', 'surface', or 'volume')    
    """
    assert plant == "maize" or plant == "soybean" or plant == "springbarley", "plant should be 'maize', or 'soybean' or 'springbarley' "
    assert dim == "3D" or dim == "1D" or dim == "2D", "dim should be '1D' or '3D'"
    assert soil in ["hydrus_loam", "hydrus_clay", "hydrus_sand", "hydrus_sandyloam"], "soil should be 'hydrus_loam', 'hydrus_clay', 'hydrus_sand' or 'hydrus_sandyloam' "
    assert outer_method in ["voronoi", "length", "surface", "volume"], "outer_method should be 'voronoi', 'length', 'surface', or 'volume'"

    # Hidden parameters
    trans_ = 0.5  # cm / day
    wilting_point = -15000  # cm

    # Numeric parameters
    random_seed = 1  # random seed
    slope = "1000"  # cm

    soil_, table_name = soil_vg_(soil)
    if plant == "maize":
        min_b, max_b, cell_number = maize_(dim)
        param_name = "Zeamays_synMRI_modified.xml"
        rs_age = 8 * 7  # 56 days
    elif plant == "soybean":
        min_b, max_b, cell_number = soybean_(dim)
        param_name = "Glycine_max_Moraes2020_opt2_modified.xml"
        rs_age = 6 * 7  # 42 days
    elif plant == "springbarley":
        min_b, max_b, cell_number = springbarley_(dim)
        param_name = "spring_barley_CF12.xml"
        rs_age = 7 * 7  # 49 days

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
    elif dim == "2D":
        s.createGrid(min_b, max_b, cell_number, periodic = True)
    elif dim == "3D":
        s.createGrid(min_b, max_b, cell_number, periodic = True)
    s.setHomogeneousIC(initial, True)  # cm pressure head top, equilibrium (contant total potential)
    s.setTopBC("noFlux")
    s.setBotBC("noFlux")
    s.setVGParameters([soil_])
    s.setParameter("Newton.EnableChop", "True")
    s.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
    s.setParameter("Soil.SourceSlope", slope)

    print("initializeProblem()")
    sys.stdout.flush()

    s.initializeProblem()
    s.setCriticalPressure(wilting_point)
    s.ddt = 1.e-5  # [day] initial Dumux time step
    print("initializeProblem() done")
    sys.stdout.flush()

    """ root hydraulic model"""
    rs = pb.MappedRootSystem()
    rs.setSeed(random_seed)
    rs.readParameters("data/" + param_name)  #

    if plant == "maize":
        params = rs.getRootRandomParameter()
        for p in params:
            p.a = 2. * p.a  # at least ..... TODO parameterisation

    if plant == "springbarley":
        print("springbarley")
        params = rs.getRootRandomParameter()
        params[2].lmax *= 2
        params[1].theta = 1.31  # why is the tap root not always 0?

    #  = rs.getRootSystemParameter()
    # seed_param.seedPos = pb.Vector3d(0., 0., -3.)  #################################################################################

    # seed = rs.getRootSystemParameter()  # SeedRandomParameter
    # seed.firstSB = 1.e6  #################################################################################

    rs.setGeometry(pb.SDF_PlantBox(1.e6, 1.e6, np.abs(min_b[2])))
    # rs.initializeDB(4, 5)
    rs.initialize()
    print("simulating root system")
    sys.stdout.flush()
    rs.simulate(rs_age, True)
    print("initializing hydraulic model")
    sys.stdout.flush()

    params = PlantHydraulicParameters()
    r = PlantHydraulicModel("Doussan", rs, params)

    r.rs.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                            pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), False)  # cutting
    if plant == "maize":
        # const_conductivities(r)
        maize_conductivities(params, 1., 1.)
    elif plant == "soybean":
        # const_conductivities(r)
        lupine_conductivities(params, 1., 1.)
    elif plant == "springbarley":
        # const_conductivities(r)
        springbarley_conductivities(params)

    print("Hydraulic model done()")

    """ coupling roots to macroscopic soil """
    if dim == "1D":
        picker = lambda x, y, z: s.pick([0., 0., z])
    elif dim == "2D":
        picker = lambda x, y, z: s.pick([x, 0., z])
    elif dim == "3D":
        picker = lambda x, y, z: s.pick([x, y, z])

    r.rs.setSoilGrid(picker)  # maps segment ############################## rs?
    sys.stdout.flush()
    seg2cell = r.rs.seg2cell
    sys.stdout.flush()
    ns = len(r.rs.segments)
    mapping = np.array([seg2cell[j] for j in range(0, ns)])
    sys.stdout.flush()

    """ outer radii """
    if outer_method == "voronoi":
        print("voronoi", outer_method)
        sys.stdout.flush()
        outer_ = PerirhizalPython(rs).get_outer_radii_bounded_voronoi()
        if np.sum([np.isnan(outer_)]) > 0:
            print("set_scenario(): NaNs in get_outer_radii_bounded_voronoi are replaced by mean", np.sum([np.isnan(outer_)]))
            outer_mean = np.nanmean(outer_) * np.ones(outer_.shape)
            outer_[np.isnan(outer_)] = outer_mean[np.isnan(outer_)]
    else:
        print("other:", outer_method)
        sys.stdout.flush()
        outer_ = PerirhizalPython(rs).get_outer_radii(outer_method)

    print("done")
    sys.stdout.flush()

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
    return RegularGridInterpolator((kx_, sx_, inner_, outer_), sra_table)  # default is 'linear' (method = 'nearest')


def write_files(file_name, hx, hsr, sink, times, trans, trans2, hs, wall_time = 0.):
    """  saves numpy arrays as npy files """
    np.save('results/hx_' + file_name, np.array(hx))  # xylem pressure head per segment [cm]
    np.save('results/hsr_' + file_name, np.array(hsr))  # pressure head at interface per segment [cm]
    np.save('results/sink_' + file_name, -np.array(sink))  # sink per soil cell [cm3/day]
    np.save('results/transpiration_' + file_name, np.vstack((times, -np.array(trans), -np.array(trans2))))  # time [day], transpiration [cm3/day]
    np.save('results/hs_' + file_name, np.array(hs))  # soil water matric potential per soil cell [cm]
    np.save('results/time_' + file_name, np.array(wall_time))

