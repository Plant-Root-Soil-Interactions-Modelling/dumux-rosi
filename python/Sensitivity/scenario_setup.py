""" 
    Functions to simplify setup of the scenarios for the INARI project
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

import os
import timeit
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); size = comm.Get_size()

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


def prepare_conductivities(mods):
    """ for conductivities from mecha, radii must be adjusted in a first step"""
    data = []
    cm = mods["conductivity_mode"]
    if cm == "from_mecha":
        ind_ = [mods["conductivity_index1"], mods["conductivity_index2"], mods["conductivity_index3"]]
        files = []
        for file in os.listdir(mods["mecha_path"]):
            if file.endswith(".npy"):  # Check if the file is a .npy file
                files.append(file)
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
    return data


def attach_conductivitis(params):
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


def write_results(file_name, pot_trans_, psi_x_, psi_i_, sink_, times_, act_trans_, psi_s_, vol_, surf_, krs_, depth_, collar_pot_):
    """  saves results from run_sra numpy arrays in a npz file 
    TODO eventually needs revisions, better naming and documentation"""
    np.savez("results/" + file_name,
             pot_trans = np.array(pot_trans_),
             psi_x = psi_x_,
             psi_rs = psi_i_,
             sink = sink_,
             times = times_,
             act_trans = act_trans_,
             psi_s = psi_s_,
             vol = vol_, surf =
             surf_, krs = krs_,
             depth = depth_,
             collar_pot = collar_pot_)


def write_cplantbox_results(file_name, length, surface, volume, depth, RLDmean, RLDz, krs, SUFz, RLD, SUF, area_, write_all = True):
    """  saves results from run_cplantbox in a npz file """
    if write_all:
        np.savez("results_cplantbox/" + file_name, length = length, surface = surface, volume = volume, depth = depth, RLDmean = RLDmean, RLDz = RLDz, krs = krs, SUFz = SUFz, area = area_, RLD = RLD, SUF = SUF)
    else:
        np.savez("results_cplantbox/" + file_name, length = length, surface = surface, volume = volume, depth = depth, RLDmean = RLDmean, RLDz = RLDz, krs = krs, SUFz = SUFz, area = area_)

