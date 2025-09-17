"""
    Dynamic:

    Creates a 4d look up table for root-soil interface potentials for specific van Genucthen parameters 
    of the envirotypes
    
    Daniel Leitner, 2025     
"""

import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

import matplotlib.pyplot as plt
import functional.van_genuchten as vg
from functional.Perirhizal import *

SMALL_SIZE = 16
MEDIUM_SIZE = 16
BIGGER_SIZE = 16
plt.rc('font', size = SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize = SMALL_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize = MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize = SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize = BIGGER_SIZE)  # fontsize of the figure title
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

if __name__ == "__main__":

    soil = {}
    soil[0] = [0.0639, 0.3698, 0.0096, 1.4646, 4.47]
    soil[1] = [0.0619, 0.3417, 0.0132, 1.3258, 2.03]
    soil[36] = [0.0760, 0.3863, 0.0091, 1.4430, 2.99]
    soil[5] = [ 0.0451, 0.3721, 0.0325, 1.4393, 34.03]
    soil[59] = [0.0534, 0.3744, 0.0171, 1.4138, 13.09]

    for pair in soil.items():
        vg.plot_retention_curve(vg.Parameters(pair[1]), label_ = pair[0])

    plt.legend()
    plt.show()

    # filename = "data/envirotype0"
    # sp = vg.Parameters(soil[0])
    # vg.create_mfp_lookup(sp)
    # peri = PerirhizalPython()
    # peri.create_lookup_mpi(filename, sp)  # takes some hours
    #
    # filename = "data/envirotype1"
    # sp = vg.Parameters(soil[1])
    # vg.create_mfp_lookup(sp)
    # peri = PerirhizalPython()
    # peri.create_lookup_mpi(filename, sp)  # takes some hours
    #
    # filename = "data/envirotype36"
    # sp = vg.Parameters(soil[36])
    # vg.create_mfp_lookup(sp)
    # peri = PerirhizalPython()
    # peri.create_lookup_mpi(filename, sp)  # takes some hours
    #
    # filename = "data/envirotype5"
    # sp = vg.Parameters(soil[5])
    # vg.create_mfp_lookup(sp)
    # peri = PerirhizalPython()
    # peri.create_lookup_mpi(filename, sp)  # takes some hours
    #
    # filename = "data/envirotype59"
    # sp = vg.Parameters(soil[59])
    # vg.create_mfp_lookup(sp)
    # peri = PerirhizalPython()
    # peri.create_lookup_mpi(filename, sp)  # takes some hours

