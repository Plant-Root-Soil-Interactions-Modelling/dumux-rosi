""" 
 tests resolution dependency of krs, of a single root 
 
 gladly, no dependency 
"""

import sys; sys.path.append("../../../modules/"); sys.path.append("../../../../../CPlantBox/");  sys.path.append("../../../../../CPlantBox/src/python_modules")
sys.path.append("../../../../build-cmake/cpp/python_binding/"); sys.path.append("../../../modules/fv/");
sys.path.append("../"); sys.path.append("../scenarios/");

import plantbox as pb  # CPlantBox
from rosi_richards import RichardsSP  # C++ part (Dumux binding), macroscopic soil model
from richards import RichardsWrapper  # Python part, macroscopic soil model
from xylem_flux import *  # root system Python hybrid solver

import vtk_plot as vp
import van_genuchten as vg
import aggregated_rs as agg

import matplotlib.pyplot as plt
import numpy as np

""" 
Parameters  
"""

""" root system """
a = 0.05  # cm

""" 
Initialize xylem model 
"""
print("Meunier")
for n in [10, 50, 100, 200]:
    r = agg.create_singleroot(ns = n, l = 100 , a = a)
    agg.init_comp_conductivities_const(r)
    krs, j = r.get_krs(0.)
    print(n, ": krs", krs, "j", j, krs / (2 * a * np.pi), krs / (2 * a * np.pi * 50))  # krs / (root surface) ~= kr
    print()

print("\nDoussan")
for n in [10, 50, 100, 200]:
    r = agg.create_singleroot(ns = n, l = 100 , a = a)
    agg.init_comp_conductivities_const(r)
    r.linearSystem = r.linearSystem_doussan  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    krs, j = r.get_krs(0.)
    print(n, ": krs", krs, "j", j, krs / (2 * a * np.pi), krs / (2 * a * np.pi * 50))  # krs / (root surface) ~= kr
    print()

    # segs = r.rs.segments
    # nodes = r.rs.nodes
    # z_ = [0.5 * (nodes[s.x].z + nodes[s.y].z) for s in segs]
