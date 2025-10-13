""" 
    Dynamic:

    Simulates dynamic water movement for a single fixed soybean scenario
    
    either reproduce a parameter set given by a json file (script top), 
    or from scratch (script bot)
    
    Daniel Leitner, 2025
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

import cProfile
import os
import json

import plantbox as pb
import visualisation.vtk_plot as vp
import functional.van_genuchten as vg

import scenario_setup as scenario
import evapotranspiration as evap
import soil_model
import hydraulic_model

import sra_new
import run_sra

import numpy as np

sim_time = 5 # 87.5 # 87.5  # 40  # 87.5  # 87.5  # 87.5  # [day]
envirotype = 0
theta1 = None  # if none leave unmodified
src = None  # if none leave unmodified

# # rerun experiment from json
# folder_path = "results_cplantbox/"
# exp_name = "soybean_all14_e6ff06381be42906c846e6a7749fbaf2da9f29d01b8a6e04a85fa32e1b1cffd1"
#
# file_path = os.path.join(folder_path, exp_name + "_mods.json")
# with open(file_path, 'r', encoding = 'utf-8') as file:
#     params = json.load(file)
# print(params)
#
# enviro_type = params["enviro_type"]
# sim_time = params["sim_time"]
# assert exp_name == params["exp_name"], "wrong file name"
# params.pop("exp_name")
# params.pop("enviro_type")
# params.pop("sim_time")
# # params["mecha_path"] = "/home/daniel/Dropbox/Code/granar/mecha_results" # local
# params["mecha_path"] = "mecha_results"  # for cluster
#
# print(exp_name, enviro_type, sim_time)
# # sim_time = 1.
#
# run_sra.run_soybean("dummy{:g}".format(envirotype), envirotype, sim_time, params, save_all = True)
#
# print("fin.")

# 1
hairsZone = 1.7
hairsLength = 0.1
hairsElongation = 0.3

# # 2
# hairsZone = 0.4
# hairsLength = 0.3
# hairsElongation = 0.3

# hairsZone = 0.
# hairsLength = 0.
# hairsElongation = 0.

# mods = { "a145": 0.2, "a2": 0.04, "a3": 0.04,
#         "hairsZone145":hairsZone,
#         "hairsZone2":hairsZone,
#         "hairsZone3":hairsZone,
#         "hairsLength145":hairsLength,
#         "hairsLength2":hairsLength,
#         "hairsLength3":hairsLength,
#         "hairsElongation": hairsElongation,
#         "dx": 0.1
#         }
#

# # 3
mods = {
        "output_times": [40, 60],

        "conductivity_mode": "scale",
        "scale_kr":1.,
        "scale_kx":1.,

        "a145": 0.2, "a2": 0.04, "a3": 0.04,
        "hairsZone145":0, "hairsZone2":1.7, "hairsZone3":0.4,
        "hairsLength145":0, "hairsLength2":0.1, "hairsLength3":0.3,
        "hairsElongation": 0.3,

        "dx": 0.1,

        "bot_bc": "noFlux",
        }



import pstats

# Creating profile object
ob = cProfile.Profile()

ob.enable()
run_sra.run_soybean("soybean_test", envirotype, sim_time, mods, save_all = True)
ob.disable()

stats = pstats.Stats(ob)
stats.sort_stats('tottime')
stats.print_stats(20)



# # kr = np.zeros((3,))
# # kr_old = np.zeros((2,))
# # kx = np.zeros((3,))
# # kx_old = np.zeros((2,))
# #
# # kr[0] = 0.1
# # kr_old[0] = 0
# # kx[0] = 0.05
# # kx_old[0] = 0.1
# #
# # mods = {
# # "filename": "data/Glycine_max_Moraes2020_singleroot.xml",
# # "initial_age": 1.,
# # "initial_totalpotential":-500
# # }
# #
# # cu = run_sra.run_soybean("soybean_test_{:g}".format(1), envirotype, sim_time, mods, kr, kx, kr_old, kx_old, save_all = True)
#
# #
# # n = len(rs.radii)
# # radii = np.array([rs.getEffectvieRadius(i) for i in range(0, n)])
# # rs.radii = radii
# #
# # ana = pb.SegmentAnalyser(rs.mappedSegments())
# # ana.addData("distanceTip", rs.distanceTip)
# # ana.addAge(sim_time)
# #
# # vp.plot_roots(ana, "distanceTip")

