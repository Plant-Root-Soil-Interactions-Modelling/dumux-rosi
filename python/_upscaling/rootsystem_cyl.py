""" 
Single root scenario - soil depletion due to sinusoidal transpiration over 21 days

coupled to cylindrical rhizosphere models using 1d axi-symetric richards equation (DUMUX solver) (cyl)
"""
import sys; sys.path.append("../modules/"); sys.path.append("../../../CPlantBox/");  sys.path.append("../../../CPlantBox/src/python_modules")
sys.path.append("../../build-cmake/cpp/python_binding/")  # DUMUX solver

from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()

import scenario_setup as scenario
from rhizo_models import RhizoMappedSegments
import cyl
import sra

""" parameters   """
min_b = [-6, -1.5, -150.]  # domain 12cm x 3cm x 150cm
max_b = [6, 1.5, 0.]
cell_number = [12, 3, 150]  # 1 cm3

theta_r = 0.025  # sandy loam
theta_s = 0.403
alpha = 0.0383  # (cm-1) soil
n = 1.3774
k_sat = 60.  # (cm d-1)
soil_ = [theta_r, theta_s, alpha, n, k_sat]

trans = 0.6 * (12 * 3)  # cm3/day

sim_time = 21  #  [day]
dt = 20 / (24 * 3600)  # time step [day]

""" rhizosphere model parameters """
wilting_point = -10000  # cm
nc = 10  # dof+1
logbase = 0.5  # according to Mai et al. (2019)
mode = "dumux"

""" initialize """
s, soil = scenario.create_soil_model(soil_, min_b, max_b, cell_number, p_top = -330, p_bot = -180)
r = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, "results/wheat.rsml")  # created by rootsystem.py
rs = RhizoMappedSegments(r, wilting_point, nc, logbase, mode)

""" sanity checks """
if rank == 0:
    r.test()  # sanity checks

""" numerical solution """
water0 = s.getWaterVolume()  # total initial water volume in domain

# psi_x_, psi_s_, sink_, x_, y_, psi_s2_ = cyl.simulate_const(s, rs, trans, sim_time, dt)
psi_x_, psi_s_, sink_, x_, y_, psi_s2_ = sra.simulate_const(s, rs, trans, sim_time, dt)

water = s.getWaterVolume()

""" output """
if rank == 0:

    scenario.write_files("rootsystem_cyl", psi_x_, psi_s_, sink_, x_, y_, psi_s2_)

    print("\ntotal uptake", water0 - water, "cm3")
    print("fin")
