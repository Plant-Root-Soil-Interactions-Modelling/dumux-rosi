""" 
Functions to simplify setup of the scenarios for the INARI project
"""

import sys;
sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../modules/");
sys.path.append("/data");
sys.path.append("../../../CPlantBox/src/python_modules");
sys.path.append("../../../CPlantBox/src/functional/");
sys.path.append("../../../CPlantBox/src/rsml/");
sys.path.append("../../../CPlantBox/src/visualisation/")
sys.path.append("../../../CPlantBox/src/structural/")
sys.path.append("../../../CPlantBox/src/external/")
sys.path.append("../../../CPlantBox/");

import numpy as np
import timeit
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
from scipy import interpolate

from rosi_richardsnc import RichardsNCSP  # C++ part (Dumux binding), macroscopic soil model
from rosi_richards import RichardsSP  # C++ part (Dumux binding), macroscopic soil model
from richards import RichardsWrapper  # Python part, macroscopic soil model

import plantbox as pb  # CPlantBox
import van_genuchten as vg
import evapotranspiration as evap
from xylem_flux import *
from datetime import *

SMALL_SIZE = 16
MEDIUM_SIZE = 16
BIGGER_SIZE = 164
plt.rc('font', size = SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize = SMALL_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize = MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize = SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize = BIGGER_SIZE)  # fontsize of the figure title
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']


def vg_SPP(i = int(1)):
    """ Van Genuchten parameter, called by maize()  """
        
    soil = {}
    soil[0] = [0.041, 0.494, 0.0256, 1.49, 245]
    soil[1] = [0.03, 0.414, 0.038, 2, 1864]
    return soil[i]


def maize_SPP(soil_= "loam"):
    """ parameters for maize simulation """
    if soil_ == "loam":
        i = 0
    else:
        i = 1
    soil = vg_SPP(i)
    min_b = [-10., -22, -74.] 
    max_b = [10., 22, 0.]
    cell_number = [10, 22, 37]
    area = 20 * 44  # cm2

    #bulk_d = {}
    #bulk_d[0] = 1.4 #g cm-3
    #bulk_d[1] = 1.5 #g cm-3
    
    Kc_value_ = {}
    Kc_value_[0] = np.array([1,1,1,1.2,1.2,1.2])
    Kc_value_[1] = np.array([1,1,1,1.07,1.2,1.2])
    Kc_value = Kc_value_[0] #2019, 2020
    Kc_days = np.array([1,42,63,98,154,288])
    
    Kc = np.zeros((Kc_days[-1]))
    dummy = 0
    for i in range(0,len(Kc)):
        if i+1 in Kc_days:
            Kc[i] = Kc_value[np.where(Kc_days==(i+1))[0]]
            dummy = dummy+1
        else:
            slope = (Kc_value[dummy]-Kc_value[dummy-1])/(Kc_days[dummy]-Kc_days[dummy-1])
            Kc[i] = Kc_value[dummy-1]+slope*((i+1)-Kc_days[dummy-1])

    #plt.plot(np.linspace(0,287,288), Kc)
            
    return soil, min_b, max_b, cell_number, area, Kc

def exudation_rates(t, comp):
    
    times = [0, 42, 63, 98, 154]
    if comp == "phenolics":
        exu_prop = [0.00012, 0.00012, 0.00003, 0.000019, 0.000015]  #[kg/(m2 day)]
    elif comp == "sugars":
        exu_prop = [0.0006, 0.0006, 0.00015, 0.00017, 0.00016]  #[kg/(m2 day)]
    elif comp == "aminoacids":
        exu_prop = [0.0001, 0.0001, 0.000025, 0.0000168, 0.0000135]  #[kg/(m2 day)]
    elif comp == "test":
        exu_prop = [0.00012, 0.00012, 0.00003, 0.000019, 0.000015]  #[kg/(m2 day)]
    else:
        print('No exudate properties found')

    f = interpolate.interp1d(times, exu_prop)  
    kex = np.array([[0., 5.], [f(t), 0.]])
        
    return kex

def init_conductivities_const(r, kr_const = 1.8e-4, kx_const = 0.1):
    """ Hydraulic conductivities  kr [1/day], kx [cm3/day] """
    r.setKr([0, kr_const, kr_const, kr_const, kr_const, kr_const])
    r.setKx([1.e3, kx_const, kx_const, kx_const, kx_const, kx_const])


def init_maize_conductivities(r, skr = 1., skx = 1.):
    """ Hydraulic conductivities for maize modified from Couvreur et al. (2012), originally from Doussan et al. (1998) """ #[age, value]

    kr00 = np.array([[0., 0.]])  # artificial shoot
    kx00 = np.array([[0., 1.e3]])  # artificial shoot

    kr0 = np.array([[-1e4, 0.], [-0.1, 0.], [0., 0.000181], [8., 0.000181], [10, 0.0000648], [18, 0.0000648], [25, 0.0000173], [300, 0.0000173]])
    kr1 = np.array([[-1e4, 0.], [-0.1, 0.], [0., 0.000181], [10., 0.000181], [16, 0.0000173], [300, 0.0000173]])

    kx0 = np.array([[0., 0.000864], [5., 0.00173], [12., 0.0295], [15., 0.0295], [20., 0.432], [300., 0.432]])
    kx1 = np.array([[0., 0.000864], [5., 0.0000864], [10., 0.0000864], [12., 0.0006048], [20., 0.0006048], [23., 0.00173], [300., 0.00173]])

    kr01 = np.minimum(skr * kr0[:, 1], 1.)
    kr11 = np.minimum(skr * kr1[:, 1], 1.)
    r.setKrTables([kr00[:, 1], kr01, kr11, kr11, kr01, kr01],
                  [kr00[:, 0], kr0[:, 0], kr1[:, 0], kr1[:, 0], kr0[:, 0], kr0[:, 0]])
    kx01 = np.minimum(skx * kx0[:, 1], 1.)
    kx11 = np.minimum(skx * kx1[:, 1], 1.)
    r.setKxTables([kx00[:, 1], kx01, kx11, kx11, kx01, kx01],
                  [kx00[:, 0], kx0[:, 0], kx1[:, 0], kx1[:, 0], kx0[:, 0], kx0[:, 0]])


def create_soil_model(soil_type, year, soil_,comp, min_b , max_b , cell_number, type, times = None, net_inf = None):
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
    s.createGrid(min_b, max_b, cell_number, False)  # [cm] #######################################################################

    # BC
    if times is not None:
        s.setTopBC("atmospheric", 0.5, [times, net_inf])  # 0.5 is dummy value
    else:
        s.setTopBC("noFlux")
    s.setBotBC("noFlux") #in acc. with Jorda et al. (2022), however, they assume inflow if h>0

    if type == 2:  # solute BC
        s.setTopBC_solute("outflow", 0.)
        s.setBotBC_solute("outflow", 0.)

# Paramters
    s.setVGParameters([soil_])
    s.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
    if type == 2:
        if comp == "phenolics": 
            s.setParameter("Component.MolarMass", "1.2e-2")  # carbon 12 g/mol
            s.setParameter("Component.LiquidDiffusionCoefficient", "2e-10")  # m2 s-1 # Srinivas et al. 2011
            if soil_type == "loam": 
                s.setParameter("Component.BufferPower", "92")  # buffer power = \rho * Kd [1]
            else:
                s.setParameter("Component.BufferPower", "44")  # buffer power = \rho * Kd [1], Cecchi et al. (2004)
            s.decay = 0.3 # day-1 Wang et al. (2016)
            
        elif comp == "sugars":
            s.setParameter("Component.MolarMass", "1.2e-2")  # carbon 12 g/mol
            s.setParameter("Component.LiquidDiffusionCoefficient", "6e-10")  # m2 s-1 # Ribeiro et al. 2006, Hartig et al. 2018
            
            s.setParameter("Component.BufferPower", "0")  # buffer power = \rho * Kd [1]
            s.decay = 16 # day-1 Gunina & Kuzyakov (2015)
            
        elif comp == "aminoacids":
            s.setParameter("Component.MolarMass", "1.2e-2")  # carbon 12 g/mol
            s.setParameter("Component.LiquidDiffusionCoefficient", "6e-10")  # m2 s-1 # Ma et al. 2005, Virk et al. 2015

            if soil_type == "loam": 
                s.setParameter("Component.BufferPower", "126")  # buffer power = \rho * Kd [1]
            else:
                s.setParameter("Component.BufferPower", "7")  # buffer power = \rho * Kd [1], Oburger et al. 2011
            s.decay = 2 # day-1 Fischer et al. (2007)
        elif comp == "test":
            s.setParameter("Component.MolarMass", "1.2e-2")  # carbon 12 g/mol
            s.setParameter("Component.LiquidDiffusionCoefficient", "6e-10")  # m2 s-1
            s.setParameter("Component.BufferPower", "0")  
            s.decay = 0 # day-1 
        else:
            print('no exudate compound defined') 
            
            
    s.initializeProblem()
    wilting_point = -15000
    s.setCriticalPressure(wilting_point)  # for boundary conditions constantFlow, constantFlowCyl, and atmospheric
    s.ddt = 1.e-5  # [day] initial Dumux time step

    # IC
    df = pd.read_csv("data_magda/init_pot_"+str(year)+".csv")  # initial potential
    h  = np.flip(df[soil_type].loc[:].values) #cm
    h = np.repeat(h[:,np.newaxis],cell_number[0],axis=1) #x-axis
    h = np.repeat(h[:,:,np.newaxis],cell_number[1],axis=2) #y-axis
    h = h.flatten()
    #h = np.ones((20*45*75))*-100 #TODO
    s.setInitialConditionHead(h)  # cm

    if type == 2:
        c = np.zeros((cell_number[0]*cell_number[1]*cell_number[2])) #TODO
        s.setInitialCondition(c, 1)  # kg/m3

    # plt.plot(h, np.linspace(-200., 0., h.shape[0]))
    # plt.xlabel("soil matric potential [cm]")
    # plt.ylabel("depth (cm)")
    # plt.tight_layout()
    # plt.show()
    # plt.plot(c, np.linspace(-200, 0., c.shape[0]))
    # plt.xlabel("nitrate concentration [g/cm3]")
    # plt.ylabel("depth (cm)")
    # plt.tight_layout()
    # plt.show()

    return s, soil

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


def create_mapped_rootsystem(min_b , max_b , cell_number, soil_model, fname, stochastic = False, mods = None):
    """ loads a rmsl file, or creates a rootsystem opening an xml parameter set,  
        and maps it to the soil_model """
    global picker  # make sure it is not garbage collected away...

    if fname.endswith(".rsml"):
        r = XylemFluxPython(fname)
    elif fname.endswith(".xml"):
        seed = 1

        rs = pb.MappedRootSystem()
        rs.setSeed(seed)
        rs.readParameters(fname)
        if not stochastic:
            set_all_sd(rs, 0.)

        rs.setGeometry(pb.SDF_PlantBox(max_b[0]*2, max_b[1]*2, np.abs(min_b[2])))
        rs.initializeLB(5, 4)
        rs.simulate(1., True)
        r = XylemFluxPython(rs)

    r.rs.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                            pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), cut = False)

    picker = lambda x, y, z: soil_model.pick([0., 0., z])  #  function that return the index of a given position in the soil grid (should work for any grid - needs testing)
    r.rs.setSoilGrid(picker)  # maps segments, maps root segements and soil grid indices to each other in both directions
    # comm.barrier()
    # print("survived setSoilGrid", rank)

    # if rank == 0:
    init_maize_conductivities(r)

    return r


def write_files(file_name, psi_x, psi_i, sink, times, trans, psi_s, vol_, surf_, krs_, depth_,  dist, con, l, conc = None, c_ = None, mass_soil_c_ = None):
    """  saves numpy arrays as npy files """

    np.save('results/psix_' + file_name, np.array(psi_x))  # xylem pressure head per segment [cm]
    np.save('results/psiinterface_' + file_name, np.array(psi_i))  # pressure head at interface per segment [cm]
    np.save('results/sink_' + file_name, -np.array(sink))  # sink per segment [cm3/day]
    np.save('results/transpiration_' + file_name, np.vstack((times, -np.array(trans))))  # time [day], transpiration [cm3/day]
    np.save('results/soil_' + file_name, np.array(psi_s))  # soil potential per cell [cm]

    np.save('results/vol_' + file_name, np.array(vol_))  # volume per subType [cm3]
    np.save('results/surf_' + file_name, np.array(surf_))  # surface per subType [cm2]
    np.save('results/krs_' + file_name, np.array(krs_))  # soil potential per cell [cm2/day]
    np.save('results/depth_' + file_name, np.array(depth_))  # root system depth [cm]

    if conc is not None:
        np.save('results/soilc_' + file_name, np.array(conc))  # solute concentration [g/cm3]
    if c_ is not None:
        np.save('results/carbon_' + file_name,  np.vstack((times, -np.array(c_))))  # exudation per segment  [g/day]

    np.savez('results/solute_' + file_name, np.array(dist), np.array(con),  np.array(l))  # distance from root [cm], solute concentrations in cylinders [g/cm3], segment length [cm]

    np.savez('results/check_solute_' + file_name, np.array(c_), np.array(mass_soil_c_))  # exudation per segment [g/d], solute mass in soil domain [g]



if __name__ == '__main__':

    pass