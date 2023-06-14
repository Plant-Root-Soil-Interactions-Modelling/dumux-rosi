import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

import plantbox as pb  # CPlantBox
from rosi_richards import RichardsSP  # C++ part (Dumux binding), macroscopic soil model
from richards import RichardsWrapper  # Python part, macroscopic soil model
from functional.xylem_flux import *  # root system Python hybrid solver
from rhizo_models import *  # Helper class for cylindrical rhizosphere models

import visualisation.vtk_plot as vp
import functional.van_genuchten as vg
from functional.root_conductivities import *

import numpy as np
import timeit
import matplotlib.pyplot as plt
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size()

from scenario_setup import *



def initialize_xylem_model(plant, wilting_point, NC, logbase, mode,min_b,max_b,cell_number,age_dependent,s):
    rs = RhizoMappedSegments(plant, wilting_point, NC, logbase, mode)
    rs.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                            pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), True)
    r = XylemFluxPython(rs)  # wrap the xylem model around the MappedSegments
    init_conductivities_growth(r, age_dependent, 0.05)  # age_dependent is a boolean, root conductivies are given in the file src/python_modules/root_conductivities.py

    picker = lambda x, y, z: s.pick([x, y, z])  #  function that return the index of a given position in the soil grid (should work for any grid - needs testing)
    rs.setSoilGrid(picker)  # maps segments, maps root segements and soil grid indices to each other in both directions
    rs.set_xylem_flux(r)
    
    # For debugging
    # r.plot_conductivities()
    r.test()  # sanity checks (todo need improvements...)
    return {'r':r,'rs':rs,'picker':picker}
