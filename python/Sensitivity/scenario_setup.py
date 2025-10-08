""" 
    Dynamic & Macroscopic:

    Functions to simplify setup of the scenarios for the INARI project
    
    e.g. set_conductivities, write_results, etc.
    
    Daniel Leitner, 2025    
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

import os
import timeit
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
# from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); size = comm.Get_size()

import plantbox as pb  # CPlantBox
import functional.van_genuchten as vg
from functional.PlantHydraulicParameters import PlantHydraulicParameters
from functional.PlantHydraulicModel import HydraulicModel_Doussan
from functional.PlantHydraulicModel import HydraulicModel_Meunier
import evapotranspiration as evap
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


def vg_enviro_type(i:int):
    """ Van Genuchten parameter for enviro-types, called by maize() and soybean() """
    soil = {}
    soil[0] = [0.0639, 0.3698, 0.0096, 1.4646, 4.47]
    soil[1] = [0.0619, 0.3417, 0.0132, 1.3258, 2.03]
    soil[36] = [0.0760, 0.3863, 0.0091, 1.4430, 2.99]
    soil[5] = [ 0.0451, 0.3721, 0.0325, 1.4393, 34.03]
    soil[59] = [0.0534, 0.3744, 0.0171, 1.4138, 13.09]
    table_name = "envirotype{:s}".format(str(i))
    return soil[i], table_name


def maize(i:int):
    """ parameters for a maize simulation for envirotye number i    
    
    soil                 list of VG parameters
    table_name           look up table for root-soil interface matric potentials 
    min_b, max_b         domain boundaries
    area                 surface area fo domain
    Kc_maize             crop coefficent (Allen, et al. 1998)
    """
    soil, table_name = vg_enviro_type(i)
    min_b = [-38., -8., -200.]  # data from INARI
    max_b = [38., 8., 0.]
    cell_number = [1, 1, 200]
    area = 76. * 16  # cm2
    Kc_maize = 1.2  # book "crop evapotranspiration"
    return soil, table_name, min_b, max_b, cell_number, area, Kc_maize


def soybean(i:int):
    """ parameters for soybean simulation for envirotype number i 
    
    soil                 list of VG parameters
    table_name           look up table for root-soil interface matric potentials 
    min_b, max_b         domain boundaries
    area                 surface area fo domain
    Kc_soybean           crop coefficent (Allen, et al. 1998)
    """
    soil, table_name = vg_enviro_type(i)
    min_b = [-38, -1.5, -200.]  # data from INARI
    max_b = [38, 1.5, 0.]
    cell_number = [1, 1, 200]
    area = 76 * 3  # cm2
    Kc_soybean = 1.15  # book "crop evapotranspiration" Allen, et al 1998
    return soil, table_name, min_b, max_b, cell_number, area, Kc_soybean


def init_conductivities_const(r, kr_const = 1.8e-4, kx_const = 0.1):
    """ Sets constant hydraulic conductivities kr [1/day] and kx [cm3/day] """
    r.set_kr_const(0., 0)  # artificial shoot
    r.set_kx_const(1.e3, 0)
    r.set_kr_const(kr_const, [1, 2, 3, 4, 5])  # subTypes = 1-5
    r.set_kx_const(kx_const, [1, 2, 3, 4, 5])


def init_lupine_conductivities3(r, ykr1, okr1, ykr2, okr2, ykr3, okr3, ykx1, okx1, ykx2, okx2, ykx3, okx3,
                                yt1, yt2, yt3, delta1 = 1, delta2 = 2, delta3 = 3):
    """ 15+(3) variables set up: young and old kr and kx (4) for 3 root types (145, 2, 3), 
        and 3 transition times (from young to old) and optional slopes """
    r.set_kr_age_dependent([0.], [0.], 0)  # artificial shoot
    r.set_kx_age_dependent([0.], [1.e3], 0)
    kx1 = np.array([[0., ykx1], [yt1, ykx1], [yt1 + delta1, okx1], [yt1 + delta1, okx1]])  # increasing with age
    kr1 = np.array([[0., ykr1], [yt1, ykr1], [yt1 + delta1, okr1], [yt1 + delta1, okr1]])  # decreasing with age
    kx2 = np.array([[0., ykx2], [yt2, ykx2], [yt2 + delta2, okx2], [yt2 + delta2, okx2]])  # increasing with age
    kr2 = np.array([[0., ykr2], [yt2, ykr2], [yt2 + delta2, okr2], [yt2 + delta2, okr2]])  # decreasing with age
    kx3 = np.array([[0., ykx3], [yt3, ykx3], [yt3 + delta3, okx3], [yt3 + delta3, okx3]])  # increasing with age
    kr3 = np.array([[0., ykr3], [yt3, ykr3], [yt3 + delta3, okr3], [yt3 + delta3, okr3]])  # decreasing with age
    r.set_kr_age_dependent(kr1[:, 0], kr1[:, 1], [1, 4, 5]),  # age, value, subType
    r.set_kx_age_dependent(kx1[:, 0], kx1[:, 1], [1, 4, 5])
    r.set_kr_age_dependent(kr2[:, 0], kr2[:, 1], 2)
    r.set_kx_age_dependent(kx2[:, 0], kx2[:, 1], 2)
    r.set_kr_age_dependent(kr3[:, 0], kr3[:, 1], 3)
    r.set_kx_age_dependent(kx3[:, 0], kx3[:, 1], 3)


# scenario.init_lupine_conductivities_sa(r, kr[0], kr_old[0], kr[1], kr_old[1], kr[2], kx[0], kx_old[0], kx[1], kx_old[1], kx[2])
def init_lupine_conductivities_sa(r, ykr1, okr1, ykr2, okr2, kr3_, ykx1, okx1, ykx2, okx2, kx3_):
    """ 10 variable set up for sensitivity analysis"""
    r.set_kr_age_dependent([0.], [0.], 0)  # artificial shoot
    r.set_kx_age_dependent([0.], [1.e3], 0)
    kx1 = np.array([[0., ykx1], [14., ykx1], [28., okx1], [28., okx1]])  # increasing with age
    kr1 = np.array([[0., ykr1], [7., ykr1], [28., okr1], [28., okr1]])  # decreasing with age
    kx2 = np.array([[0., ykx2], [7., ykx2], [14., okx2], [14., okx2]])  # increasing with age
    kr2 = np.array([[0., ykr2], [7., ykr2], [14., okr2], [14., okr2]])  # decreasing with age
    kr3 = np.array([[0., kr3_], [100., kr3_]])  # constant
    kx3 = np.array([[0., kx3_], [100., kx3_]])  # constant
    r.set_kr_age_dependent(kr1[:, 0], kr1[:, 1], [1, 4, 5]),  # age, value, subType
    r.set_kx_age_dependent(kx1[:, 0], kx1[:, 1], [1, 4, 5])
    r.set_kr_age_dependent(kr2[:, 0], kr2[:, 1], 2)
    r.set_kx_age_dependent(kx2[:, 0], kx2[:, 1], 2)
    r.set_kr_age_dependent(kr3[:, 0], kr3[:, 1], 3)
    r.set_kx_age_dependent(kx3[:, 0], kx3[:, 1], 3)


def init_lupine_conductivities(r, skr = 1., skx = 1.):
    """ 2 variables set up: Hydraulic conductivities for lupine following Zarebanadkouki et al. (2016), two variables to scale the inputs """
    r.set_kr_age_dependent([0.], [0.], 0)  # artificial shoot
    r.set_kx_age_dependent([0.], [1.e3], 0)
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
    kx01 = np.minimum(skx * kx0[:, 1], 1.)
    kx11 = np.minimum(skx * kx1[:, 1], 1.)
    r.set_kr_age_dependent(kr0[:, 0], kr01, [1, 4, 5]),  # age, value, subType
    r.set_kx_age_dependent(kx0[:, 0], kx01, [1, 4, 5])
    r.set_kr_age_dependent(kr1[:, 0], kr11, [2, 3])
    r.set_kx_age_dependent(kx1[:, 0], kx11, [2, 3])


def sort_indices(mods):
    """ converts the file number index to anatomy index and 
    sorts according to radius to know which index applies to which subType"""
    data, ana_ind = [], []
    cm = mods["conductivity_mode"]
    if cm == "from_mecha":
        ind_ = np.array([mods["conductivity_index1"], mods["conductivity_index2"], mods["conductivity_index3"]])
        files = []
        for file in os.listdir(mods["mecha_path"]):
            if file.endswith(".npy"):  # Check if the file is a .npy file
                files.append(file)
        for i in ind_:
            file = files[int(i)]
            file_path = os.path.join(mods["mecha_path"], file)
            data.append(np.load(file_path))
            ana_ind.append(file.split("shiny")[1].split(".")[0])
    
        data = np.array(data)
        a_ = data[:, 2, 2]  # [cm] (see granar/functions.py)
        a = np.array([float(x) for x in a_])
        I = np.argsort(-a)  # descending regarding radius a
        ana_ind = np.array(ana_ind)
        return ana_ind[I]



def prepare_conductivities(mods):
    """ for conductivities from mecha, radii must be adjusted in a first step"""
    data = []
    cm = mods["conductivity_mode"]
    if cm == "from_mecha":
        ind_ = [mods["conductivity_index1"], mods["conductivity_index2"], mods["conductivity_index3"]]
        print("indices are", ind_)
        files = []
        for file in os.listdir(mods["mecha_path"]):
            if file.endswith(".npy"):  # Check if the file is a .npy file
                files.append(file)
        print(files)
        for i in ind_:
            file = files[int(i)]
            file_path = os.path.join(mods["mecha_path"], file)
            data.append(np.load(file_path))
        data = np.array(data)
        a_ = data[:, 2, 2]  # [cm] (see granar/functions.py)
        a = np.array([float(x) for x in a_])
        I = np.argsort(-a)  # descending regarding radius a
        a = a[I]
        mods["a145_a"] = a[0]
        mods["a2_a"] = a[1]
        mods["a3_a"] = a[2]
    return data


def attach_conductivitis(mods):
    """ for conductivities from mecha, radii must be adjusted in a first step"""
    assert mods["conductivity_mode"] == "from_mecha", "scenario_setup.attach_conductivitis() conductivity_mode must be 'from_mecha'"
    ind_ = [mods["conductivity_index1"], mods["conductivity_index2"], mods["conductivity_index3"]]
    files = []
    for file in os.listdir(mods["mecha_path"]):
        if file.endswith(".npy"):  # Check if the file is a .npy file
            files.append(file)
    data = []
    for i in ind_:
        file = files[int(i)]
        # ind = file.split("shiny")[1].split(".")[0]
        file_path = os.path.join(mods["mecha_path"], file)
        data.append(np.load(file_path))
    data = np.array(data)
    a_ = data[:, 2, 2]  # [cm] (see granar/functions.py)
    a = np.array([float(x) for x in a_])
    I = np.argsort(-a)  # descending regarding radius a
    a = a[I]
    mods["a145_a"] = a[0]
    mods["a2_a"] = a[1]
    mods["a3_a"] = a[2]
    ykx_ = data[:, 0, 0]  # maturity level 0
    ykr_ = data[:, 0, 1]
    okx_ = data[:, 2, 0]  # maturity level 2
    okr_ = data[:, 2, 1]
    okx = np.array([float(x) for x in okx_])  #  [cm^4/hPa/d] ~ cm3 /day (see granar/functions.py)
    okr = np.array([float(x) for x in okr_])  # [cm/hPa/d] ~ 1 /day (see granar/functions.py)
    ykx = np.array([float(x) for x in ykx_])
    ykr = np.array([float(x) for x in ykr_])
    mods["kr_young145"] = ykr[0]
    mods["kr_young2"] = ykr[1]
    mods["kr_young3"] = ykr[2]
    mods["kr_old145"] = okr[0]
    mods["kr_old2"] = okr[1]
    mods["kr_old3"] = okr[2]
    mods["kx_young145"] = ykx[0]
    mods["kx_young2"] = ykx[1]
    mods["kx_young3"] = ykx[2]
    mods["kx_old145"] = okx[0]
    mods["kx_old2"] = okx[1]
    mods["kx_old3"] = okx[2]
    return data


def set_conductivities(params, mods, data):
    """ set the hydraulic conductivities @params given in the model parameters in @param mods """

    cm = mods["conductivity_mode"]
    if cm == "age_dependent":

        init_lupine_conductivities_sa(params, mods["kr_young1"], mods["kr_old1"], mods["kr_young2"], mods["kr_old2"], mods["kr3"],
                                               mods["kx_young1"], mods["kx_old1"], mods["kx_young2"], mods["kx_old2"], mods["kx3"])
        mods.pop("kr_young1")
        mods.pop("kr_old1")
        mods.pop("kr_young2")
        mods.pop("kr_old2")
        mods.pop("kr3")
        mods.pop("kx_young1")
        mods.pop("kx_old1")
        mods.pop("kx_young2")
        mods.pop("kx_old2")
        mods.pop("kx3")

    elif cm == "from_mecha":

        yt1, yt2, yt3 = mods["conductivity_age1"], mods["conductivity_age2"], mods["conductivity_age3"]
        ykx_ = data[:, 0, 0]  # maturity level 0
        ykr_ = data[:, 0, 1]
        ya_ = data[:, 0, 2]
        kx_ = data[:, 2, 0]  # maturity level 2
        kr_ = data[:, 2, 1]
        a_ = data[:, 2, 2]
        kx = np.array([float(x) for x in kx_])  #  [cm^4/hPa/d] ~ cm3 /day (see granar/functions.py)
        kr = np.array([float(x) for x in kr_])  # [cm/hPa/d] ~ 1 /day (see granar/functions.py)
        a = np.array([float(x) for x in a_])
        ykx = np.array([float(x) for x in ykx_])
        ykr = np.array([float(x) for x in ykr_])
        ya = np.array([float(x) for x in ya_])
        I = np.argsort(-a)  # descending regarding radius a
        init_lupine_conductivities3(params, ykr[0], kr[0], ykr[1], kr[1], ykr[2], kr[2],
                                            ykx[0], kx[0], ykx[1], kx[1], ykx[2], kx[2],
                                            yt1, yt2, yt3)
        mods.pop("mecha_path")
        mods.pop("conductivity_index1")
        mods.pop("conductivity_index2")
        mods.pop("conductivity_index3")
        mods.pop("conductivity_age1")
        mods.pop("conductivity_age2")
        mods.pop("conductivity_age3")

    elif cm == "scale":

        init_lupine_conductivities(params, mods["scale_kr"], mods["scale_kx"])
        mods.pop("scale_kr")
        mods.pop("scale_kx")

    else:
        raise "run_cplantbox.run_soybean() conductivity_mode unknown"
    mods.pop("conductivity_mode")

def write_sa_numpy(file_name, sa):
    """ saves a segment analyser using numpy arrays """
    segs_ = sa.segments        
    segs = np.array([[s.x, s.y] for s in segs_], dtype = np.int64)
    nodes_ = sa.nodes
    nodes = np.array([[n.x, n.y, n.z] for n in nodes_])
    data = sa.data
    data["segments"] = segs
    data["nodes"] = nodes
    np.savez_compressed(file_name, **data)
    
def open_sa_numpy(file_name):
    """ opens a segment analyser (see write_sa_numpy) """    
    data = np.load(file_name)
    segs_ = data["segments"]
    nodes_ = data["nodes"]
    segCTs = data["creationTime"]
    radii = data["radius"]
    segs = [pb.Vector2i(s[0], s[1]) for s in segs_]
    nodes = [pb.Vector3d(n[0],n[1],n[2]) for n in nodes_]
    ana = pb.SegmentAnalyser()
    ana.segments = segs # , nodes, segCTs, radii)        
    ana.nodes = nodes
    for key, value in data.items():
        if not key in ['segments', 'nodes']:
            # print("setting", key, len(value))        
            ana.addData(key, value)    
    return ana

def write_results(file_name, r1, r2, r3 = None):
    """  saves results from run_sra numpy arrays in a npz file 
    List 1 (r1)
    times_                       times (simulation time plus initial age)
    pot_trans_                   potential transpiration [cm3/day]
    act_trans_                   actual transpiration [cm3/day]
    collar_pot_                  root collar potential [cm]  

    List 2 (r2) for each 10th time step
    times_lr_                    times for the following arrays  
    sink_                        water uptake (per cell)  
    psi_s_                       soil matric potentials (per cell)   
    net_change_                  net domain water change (including plant uptake)
    vol_                         root system volume [cm3] per subType
    surf_                        root system surface [cm2] per subType
    depth_                       root system depth [cm] 
    krs_                         root system hydraulic conductivity [cm2/day]

    List 3 (r3) of defined output times
    out_times_                   output times including start and final simulation time (plus initial age)
    ana_                         list of segment analysers at the out_times      
   """
    if r3:
        times_sa = r3[0]
    else:
        times_sa = []
   
    np.savez_compressed("results/" + file_name,
             times = r1[0], pot_trans = r1[1], act_trans = r1[2], collar_pot = r1[3],
             times_lr = r2[0], sink = r2[1], psi_s = r2[2], net_change = r2[3], 
             vol = r2[4], surf = r2[5], depth = r2[6], krs = r2[7], times_sa = times_sa)
    if r3: 
        for i,ana in enumerate(r3[1]):
           write_sa_numpy("results/sa_" + file_name +"_"+str(i), ana) 


def write_cplantbox_results(file_name, length, surface, volume, depth, RLDmean, RLDz, krs, SUFz, RLD, SUF, area_, carbon, write_all = True):
    """  saves results from run_cplantbox in a npz file """
    if write_all:
        np.savez("results_cplantbox/" + file_name, length = length, surface = surface, volume = volume, depth = depth,
                 RLDmean = RLDmean, RLDz = RLDz, krs = krs, SUFz = SUFz, area = area_, carbon = carbon, RLD = RLD, SUF = SUF)
    else:
        np.savez("results_cplantbox/" + file_name, length = length, surface = surface, volume = volume, depth = depth,
                 RLDmean = RLDmean, RLDz = RLDz, krs = krs, SUFz = SUFz, carbon = carbon, area = area_)


def get_anatomy(index):
    """ retrieves root anatomical information from the simulation index """
    # Parameter space of the Granar & Mecha simukation
    # Xylem
    xylem_vessel_diameter = np.linspace(0.02, 0.1, num = 5)  # seq(0.02, 0.1, length.out=5)
    xylem_vessel_diameter = xylem_vessel_diameter[1:4]  # skip first and last (Python slicing excludes last index)
    xylem_number_of_poles = [4]
    # Stele
    stele_cell_diameter = [0.004]  # like c(0.004)
    stele_diameter = np.linspace(0.07, 0.042, num = 5)  # seq(0.07, 0.042, length.out=5)
    stele_diameter = stele_diameter[1:4]
    # Cortex
    cortex_cell_diameter = np.linspace(0.01, 0.06, num = 5)
    cortex_cell_diameter = cortex_cell_diameter[1:4]
    cortex_layers = np.linspace(1, 15, num = 8)
    cortex_layers = cortex_layers[1:8]  # R 2:8 → Python slice 1:8
    aerenchyma_number_of_areas = [1, 4, 8, 12]  # Aerenchyma
    aerenchyma_percentage = [0., 0.07, 0.15, 0.22, 0.30]

      # c = 0
      # for (ap in aerenchyma_percentage) {
      #     for (anoa in aerenchyma_number_of_areas) {
      #       for (cl in cortex_layers) {
      #         for (ccd in cortex_cell_diameter) {
      #           for (sd in stele_diameter) {
      #             for (scd in stele_cell_diameter) {
      #               for (xnop in xylem_number_of_poles) {
      #                 for (xvd in xylem_vessel_diameter) {
      #                    c++
    index = int(index)
    n = len(aerenchyma_number_of_areas) * len(cortex_layers) * len(cortex_cell_diameter) * len(stele_diameter) * len(stele_cell_diameter) * len(xylem_number_of_poles) * len(xylem_vessel_diameter)

    aep_ind = index // n
    index = index % n
    n = len(cortex_layers) * len(cortex_cell_diameter) * len(stele_diameter) * len(stele_cell_diameter) * len(xylem_number_of_poles) * len(xylem_vessel_diameter)
    ana_ind = index // n
    index = index % n
    n = len(cortex_cell_diameter) * len(stele_diameter) * len(stele_cell_diameter) * len(xylem_number_of_poles) * len(xylem_vessel_diameter)
    cl_ind = index // n
    index = index % n
    n = len(stele_diameter) * len(stele_cell_diameter) * len(xylem_number_of_poles) * len(xylem_vessel_diameter)
    ccl_ind = index // n
    index = index % n
    n = len(stele_cell_diameter) * len(xylem_number_of_poles) * len(xylem_vessel_diameter)
    sd_ind = index // n
    index = index % n
    n = len(xylem_number_of_poles) * len(xylem_vessel_diameter)
    scd_ind = index // n
    index = index % n
    n = len(xylem_vessel_diameter)
    xnp_ind = index // n
    index = index % n
    xvd_ind = index

    r = np.array([ aerenchyma_percentage[aep_ind], aerenchyma_number_of_areas[ana_ind],
         cortex_layers[cl_ind], cortex_cell_diameter[ccl_ind],
         stele_diameter[sd_ind], stele_cell_diameter[scd_ind],
         xylem_number_of_poles[xnp_ind], xylem_vessel_diameter[xvd_ind]])

    return r



            
        
        
if __name__ == "__main__":
    print(get_anatomy(3000))
