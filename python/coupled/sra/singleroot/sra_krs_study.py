""" 
 tests resolution dependency of krs, of a single root 
 
 gladly, no dependency 
"""

import sys; sys.path.append("../../../modules/"); sys.path.append("../../../../../CPlantBox/");  sys.path.append("../../../../../CPlantBox/src/python_modules")
sys.path.append("../../../../build-cmake/cpp/python_binding/"); sys.path.append("../../../modules/fv/");
sys.path.append("../");

import plantbox as pb  # CPlantBox
from rosi_richards import RichardsSP  # C++ part (Dumux binding), macroscopic soil model
from richards import RichardsWrapper  # Python part, macroscopic soil model
from xylem_flux import *  # root system Python hybrid solver
from rhizo_models import *  # Helper class for cylindrical rhizosphere models

import vtk_plot as vp
import van_genuchten as vg
from root_conductivities import *

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

for n in [10, 50, 100, 200]:

    radii = np.array([a] * n)
    nodes = [pb.Vector3d(0, 0, 0)]
    segs = []
    for i in range(0, n):
        nodes.append(pb.Vector3d(0, 0, -(i / (n - 1)) * 50 - 0.5))
        segs.append(pb.Vector2i(i, i + 1))

    ms = pb.MappedSegments(nodes, segs, radii)
    r = XylemFluxPython(ms)  # wrap the xylem model around the MappedSegments
    init_singleroot_contkrkx(r)
    if n == 10:  # first run
        print("kr: ", r.kr_f(0., 0))
        print("kx: ", r.kx_f(0., 0))

    # r.test()  # sanity checks
    krs, j = r.get_krs(0.)
    print(n, ": krs", krs, "j", j, krs / (2 * a * np.pi), krs / (2 * a * np.pi * 50))  # krs / (root surface) ~= kr
    print()

