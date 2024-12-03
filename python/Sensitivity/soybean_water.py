""" 
    Simulates water movement for a single fixed soybean scenario
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

sim_time = 3  # 87.5  # 87.5  # 87.5  # [day]
envirotype = 0
theta1 = None  # if none leave unmodified
src = None  # if none leave unmodified

a = 0.3  # cm
hairsZone = 1.7
hairsLength = 0.1
hairsElongation = 0.3

mods = {"a145": 0.2, "a2": 0.04, "a3": 0.04,
        "hairsZone145":hairsZone, "hairsZone2":hairsZone, "hairsZone3":hairsZone,
        "hairsLength145":hairsLength, "hairsLength2":hairsLength, "hairsLength3":hairsLength,
        "hairsElongation": hairsElongation, "dx": 0.1 }

s, r = run_sra.run_soybean("soybean_test_{:g}".format(envirotype), envirotype, sim_time, mods, 1., 1., save_all = True)

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

