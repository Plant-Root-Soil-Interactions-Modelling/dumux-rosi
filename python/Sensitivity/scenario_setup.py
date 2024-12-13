""" 
Functions to simplify setup of the scenarios for the INARI project
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); size = comm.Get_size()
import numpy as np
import timeit
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

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
    """ parameters for soybean simulation for envirotye number i 
    
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
    r.setKr([0, kr_const, kr_const, kr_const, kr_const, kr_const])
    r.setKx([1.e3, kx_const, kx_const, kx_const, kx_const, kx_const])


def init_conductivities_const_growth(r, kr_const = 1.8e-4, kx_const = 0.1):
    """ Sets constant hydraulic conductivities kr [1/day] and kx [cm3/day] 
        with kr = 0 for negative ages (to mimic root growth) 
    """
    kr = np.array([[-1e4, 0.], [-0.1, 0.], [0., kr_const], [1.e4, kr_const]])
    kx = np.array([[0, kx_const], [1e4, kx_const]])
    kr00 = np.array([[0., 0.]])  # artificial shoot
    kx00 = np.array([[0., 1.e3]])  # artificial shoot
    r.setKrTables([kr00[:, 1], kr[:, 1], kr[:, 1], kr[:, 1], kr[:, 1], kr[:, 1]],
                  [kr00[:, 0], kr[:, 0], kr[:, 0], kr[:, 0], kr[:, 0], kr[:, 0]])  # for each subtype
    r.setKxTables([kx00[:, 1], kx[:, 1], kx[:, 1], kx[:, 1], kx[:, 1], kx[:, 1]],
                  [kx00[:, 0], kx[:, 0], kx[:, 0], kx[:, 0], kx[:, 0], kx[:, 0]])  # for each subtype


def init_dynamic_simple_growth(r, kr0, kr1, kx0, kx1, dt0 = 14., dt1 = 7., kr_f = 0.25, kx_f = 5.):
    """ Sets constant hydraulic conductivities kr [1/day] and kx [cm3/day] 
        with kr = 0 for negative ages (to mimic root growth) 
        with an artificial shoot for subTYpe 0
        """
    kr00 = np.array([[0., 0.]])  # artificial shoot
    kx00 = np.array([[0., 1.e3]])  # artificial shoot
    kr0 = np.array([[-1e4, 0.], [-0.1, 0.], [0., kr0], [dt0, kr_f * kr0]])  # primals
    kr1 = np.array([[-1e4, 0.], [-0.1, 0.], [0., kr1], [dt1, kr_f * kr1]])  # laterals
    kx0 = np.array([[0., kx0], [dt0, kx_f * kx0]])  # primals
    kx1 = np.array([[0., kx1], [dt1, kx_f * kx1]])  # laterals
    r.setKrTables([kr00[:, 1], kr0[:, 1], kr1[:, 1], kr1[:, 1], kr0[:, 1], kr0[:, 1]],
                  [kr00[:, 0], kr0[:, 0], kr1[:, 0], kr1[:, 0], kr0[:, 0], kr0[:, 0]])
    r.setKxTables([kx00[:, 1], kx0[:, 1], kx1[:, 1], kx1[:, 1], kx0[:, 1], kx0[:, 1]],
                  [kx00[:, 0], kx0[:, 0], kx1[:, 0], kx1[:, 0], kx0[:, 0], kx0[:, 0]])

# def init_timing(r, kr0 = 1.e-2, kx0 = 1.e-1, dt = 1.):


def init_maize_conductivities(r, skr = 1., skx = 1.):
    """ Hydraulic conductivities for maize modified from Couvreur et al. (2012), originally from Doussan et al. (1998) """

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


def init_maize_conductivities2(r, skr = 1., skx = 1.):
    """ Hydraulic conductivities for maize following Meunier et al (2018) """

    kr00 = np.array([[0., 0.]])  # artificial shoot
    kx00 = np.array([[0., 1.e3]])  # artificial shoot

    kr0 = np.array([[-1e4, 0.], [-0.1, 0.], [0., 0.47e-4], [19., 0.47e-4], [19.5, 0.43e-4], [30., 0.28e-4], [300., 0.28e-4]])  # seminal
    kr5 = np.array([[-1e4, 0.], [-0.1, 0.], [0., 0.39e-4], [13., 0.25e-4], [13.1, 0.17e-4], [30., 0.17e-4], [300., 0.17e-4]])  # crown
    # kr_brace = np.array([[-1e4, 0.], [-0.1, 0.], [0., 0.48e-4], [19., 0.48e-4], [20., 0.028e-4], [30., 0.02e-4], [300., 0.02e-4]])  # brace
    kr12 = np.array([[-1e4, 0.], [-0.1, 0.], [0., 5.62e-5], [1.e4, 5.62e-5]])  # lateral

    kx0 = np.array([[-1e4, 0.], [-0.1, 0.], [0., 1.76e-4], [7., 1.91e-4], [10., 5.45e-4], [21., 14.85e-4], [40., 15.85e-4], [300., 15.85e-4]])  # seminal
    kx5 = np.array([[-1e4, 0.], [-0.1, 0.], [0., 0.82e-4], [11.2, 3.71e-4], [28.7, 24.44e-4], [34.3, 73.63e-4], [40., 73.63e-4], [300., 73.63e-4]])  # crown
    # kx_brace = np.array([[-1e4, 0.], [-0.1, 0.], [0., 0.78e-4], [16.9, 2.56e-4], [27., 121.8e-4], [32.6, 195.6e-4], [40., 225.3e-4], [300., 225.3e-4]]) # brace
    kx12 = np.array([[0., 0.62e-4], [5.5, 1.44e-4], [7.92, 4.22e-4], [15., 4.95e-4], [300., 4.95e-4]])  # lateral

    r.setKrTables([kr00[:, 1], skr * kr0[:, 1], skr * kr12[:, 1], skr * kr12[:, 1], skr * kr0[:, 1], skr * kr5[:, 1]],
                  [kr00[:, 0], kr0[:, 0], kr12[:, 0], kr12[:, 0], kr0[:, 0], kr5[:, 0]])
    r.setKxTables([kx00[:, 1], skx * kx0[:, 1], skx * kx12[:, 1], skx * kx12[:, 1], skx * kx0[:, 1], skx * kx5[:, 1]],
                  [kx00[:, 0], kx0[:, 0], kx12[:, 0], kx12[:, 0], kx0[:, 0], kx5[:, 0]])


def init_lupine_conductivities(r, skr = 1., skx = 1.):
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
    r.params.setKrTables([kr00[:, 1], kr01, kr11, kr11, kr01, kr01],
                  [kr00[:, 0], kr0[:, 0], kr1[:, 0], kr1[:, 0], kr0[:, 0], kr0[:, 0]])
    # r.params.setKrTables([kr01, kr11, kr11, kr01, kr01],
    #               [kr0[:, 0], kr1[:, 0], kr1[:, 0], kr0[:, 0], kr0[:, 0]])

    kx01 = np.minimum(skx * kx0[:, 1], 1.)
    kx11 = np.minimum(skx * kx1[:, 1], 1.)
    r.params.setKxTables([kx00[:, 1], kx01, kx11, kx11, kx01, kx01],
                  [kx00[:, 0], kx0[:, 0], kx1[:, 0], kx1[:, 0], kx0[:, 0], kx0[:, 0]])
    # r.params.setKxTables([kx01, kx11, kx11, kx01, kx01],
    #               [kx0[:, 0], kx1[:, 0], kx1[:, 0], kx0[:, 0], kx0[:, 0]])


# scenario.init_lupine_conductivities_sa(r, kr[0], kr_old[0], kr[1], kr_old[1], kr[2], kx[0], kx_old[0], kx[1], kx_old[1], kx[2])
def init_lupine_conductivities_sa(r, ykr1, okr1, ykr2, okr2, kr3_, ykx1, okx1, ykx2, okx2, kx3_):
    """ 10 variable set up for sensitivity analysis"""
    kr00 = np.array([[0., 0.]])  # artificial shoot
    kx00 = np.array([[0., 1.e3]])  # artificial shoot

    kx1 = np.array([[0., ykx1], [14., ykx1], [28., okx1], [28., okx1]])  # increasing with age
    kx2 = np.array([[0., ykx2], [7., ykx2], [14., okx2], [14., okx2]])  # increasing with age

    kr1 = np.array([[0., ykr1], [7., ykr1], [28., okr1], [28., okr1]])  # decreasing with age
    kr2 = np.array([[0., ykr2], [7., ykr2], [14., okr2], [14., okr2]])  # decreasing with age

    kr3 = np.array([[0., kr3_], [100., kr3_]])  # constant
    kx3 = np.array([[0., kx3_], [100., kx3_]])  # constant

    r.params.setKrTables([kr00[:, 1], kr1[:, 1], kr2[:, 1], kr3[:, 1], kr1[:, 1], kr1[:, 1]],
                         [kr00[:, 0], kr1[:, 0], kr2[:, 0], kr3[:, 0], kr1[:, 0], kr1[:, 0]])

    r.params.setKxTables([kx00[:, 1], kx1[:, 1], kx2[:, 1], kx3[:, 1], kx1[:, 1], kx1[:, 1]],
                         [kx00[:, 0], kx1[:, 0], kx2[:, 0], kx3[:, 0], kx1[:, 0], kx1[:, 0]])


def init_lupine_conductivities2(r, skr = 1., skx = 1.):
    """ Hydraulic conductivities for maize following Meunier et al. (2018) !altered for tap and seminal! """

    kr00 = np.array([[0., 0.]])  # artificial shoot
    kx00 = np.array([[0., 1.e3]])  # artificial shoot

    kr = np.array([[-1e4, 0.], [-0.1, 0.], [0., 4.32e-4], [1.e4, 4.32e-4]])  # lateral
    kx = np.array([[-1e4, 0.], [-0.1, 0.], [0., 0.62e-4], [5.5, 1.44e-4], [7.92, 4.22e-4], [15., 4.95e-4], [300., 4.95e-4]])  # lateral

    r.setKrTables([kr00[:, 1], 0.1 * skr * kr[:, 1], skr * kr[:, 1], skr * kr[:, 1], 0.1 * skr * kr[:, 1], 0.1 * skr * kr[:, 1]],
                  [kr00[:, 0], kr[:, 0], kr[:, 0], kr[:, 0], kr[:, 0], kr[:, 0]])
    r.setKxTables([kx00[:, 1], 20 * skx * kx[:, 1], skx * kx[:, 1], skx * kx[:, 1], 10 * skx * kx[:, 1], 10 * skx * kx[:, 1]],
                  [kx00[:, 0], kx[:, 0], kx[:, 0], kx[:, 0], kx[:, 0], kx[:, 0]])


def write_results(file_name, pot_trans_, psi_x_, psi_i_, sink_, times_, act_trans_, psi_s_, vol_, surf_, krs_, depth_, collar_pot_):
    """  saves numpy arrays in a npz file """
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

