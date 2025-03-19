""" 
    Simulates water movement for a single fixed soybean scenario [TODO update to new conductivity parameters]
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

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

sim_time = 10  # 87.
enviro_type = 0
file_name = "singleroot_test"

kr = np.zeros((3,))
kr_old = np.zeros((2,))
kx = np.zeros((3,))
kx_old = np.zeros((2,))

# # dict = {'target': np.float64(-82.46445584745524), 'params': {'kr_old': np.float64(0.5116869897382039), 'kr_young': np.float64(0.0034975618797799743), 'kx_old': np.float64(0.023425562579446186), 'kx_young': np.float64(0.9137659997167616)}}
# dict_ = {'target': np.float64(-2389.686317347679), 'params': {'kr_old': np.float64(0.8738806656367334), 'kr_young': np.float64(0.9174124387499213), 'kx_old': np.float64(0.13290837152096602), 'kx_young': np.float64(0.9973715609813951)}}  # 80 days
# dict_ = {'target': np.float64(-1158.4424888284882), 'params': {'kr_old': np.float64(0.679069157496593), 'kr_young': np.float64(0.9186018593757479), 'kx_old': np.float64(0.0004030244893325418), 'kx_young': np.float64(0.9767591722719107)}}
dict_ = {'target': np.float64(41.86134877826826), 'params': {'kr_old': np.float64(0.9501761690709605), 'kr_young': np.float64(0.5566536315419682), 'kx_old': np.float64(0.9156064341599248), 'kx_young': np.float64(0.6415665673801286)}}
p = dict_["params"]
kr[0] = p["kr_young"] * 0.1
kr_old[0] = p["kr_old"] * 0.1
kx[0] = p["kx_young"]
kx_old[0] = p["kx_old"]

# kx = np.array([0.1, 1.e-3, 1.e-3])
# kx_old = np.array([0.35, 0.015])
# kr = np.array([1.e-3, 4.e-3, 4.e-3])
# kr_old = np.array([5e-4, 0.0015])

mods = {
"filename": "data/Glycine_max_Moraes2020_singleroot.xml",
"initial_age": 1.,
"initial_totalpotential":-100,
"domain_size": [1., 1., 200.],
"bot_bc": "noFlux"
}

cu = run_sra.run_soybean(file_name, enviro_type, sim_time, mods, kr, kx, kr_old, kx_old, save_all = True)

print("cumulative utpake:", cu)

# sim_time = 11.5  # 87.5  # 87.5  # 87.5  # [day]
# envirotype = 0
# theta1 = None  # if none leave unmodified
# src = None  # if none leave unmodified
# hairsZone = 1.
# hairsLength = 0.2
# hairsElongation = 0.3
#
# # hairsZone = 0.
# # hairsLength = 0.
# # hairsElongation = 0.
#
# mods = { "filename": "data/Glycine_max_Moraes2020_singleroot.xml",
#         "a145": 0.2, "a2": 0.04, "a3": 0.04,
#         "hairsZone145":hairsZone, "hairsZone2":hairsZone, "hairsZone3":hairsZone,
#         "hairsLength145":hairsLength, "hairsLength2":hairsLength, "hairsLength3":hairsLength,
#         "hairsElongation": hairsElongation,
#         "dx": 0.1,
#         "domain_size": [1., 1., 200.] }
#
# # mods = {"filename": "data/Glycine_max_Moraes2020_singleroot.xml",
# #         "initial_age": 1.,
# #         "initial_totalpotential":-100}
#
# s, r = run_sra.run_soybean("soybean_testsingle_{:g}".format(envirotype), envirotype, sim_time, mods, 1., 1., save_all = True)

# rs = r.ms
#
# n = len(rs.radii)
# radii = np.array([rs.getEffectvieRadius(i) for i in range(0, n)])
# rs.radii = radii
#
# ana = pb.SegmentAnalyser(rs.mappedSegments())
# ana.addData("distanceTip", rs.distanceTip)
# ana.addAge(sim_time)
#
# vp.plot_roots(ana, "distanceTip")

