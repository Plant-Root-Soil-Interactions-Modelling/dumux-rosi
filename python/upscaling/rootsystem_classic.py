""" 
Classic coupling (wheat root system) 

    * Root system potentials by meunier et al.
    * Classic coupling: root soil interface equals potential of finite volume cell
    
    works with MPI
"""
import sys; sys.path.append("../modules/"); sys.path.append("../../../CPlantBox/");  sys.path.append("../../../CPlantBox/src/python_modules")
sys.path.append("../../build-cmake/cpp/python_binding/")  # DUMUX solver

import numpy as np
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()

import scenario_setup as scenario
import vtk_plot as vp

""" parameters """
min_b = [-6, -1.5, -150.]  # domain 12cm x 3cm x 150cm
max_b = [6, 1.5, 0.]
cell_number = [12, 3, 150]  # 1 cm3

theta_r = 0.025
theta_s = 0.403
alpha = 0.0383  # (cm-1) soil
n = 1.3774
k_sat = 60.  # (cm d-1)
soil_ = [theta_r, theta_s, alpha, n, k_sat]

trans = 0.6 * (12 * 3)  # cm3/day

sim_time = 21  #  [day]
dt = 60. / (24 * 3600)  # time step [day]

""" Initialize macroscopic soil model """
s, soil = scenario.create_soil_model(soil_, min_b, max_b, cell_number, p_top = -330, p_bot = -180)
r = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, "results/wheat.rsml")  # created by rootsystem.py

""" sanity checks """
if rank == 0:
    pass
    # r.test()

""" numerical solution """
water0 = s.getWaterVolume()  # total initial water volume in domain

psi_x_, psi_s_, sink_, x_, y_, psi_s2_ = scenario.simulate_const(s, r, trans, sim_time, dt)

water = s.getWaterVolume()

""" output """
if rank == 0:

    scenario.write_files("rootsystem_classic", psi_x_, psi_s_, sink_, x_, y_, psi_s2_)

    print("\ntotal uptake", water0 - water, "cm3")
    print("fin")

    # """ VTK visualisation (does not work for MPI)"""
    # periodic = True
    # vp.plot_roots_and_soil(r.rs, "pressure head", psi_x_[-1], s, periodic, np.array(min_b), np.array(max_b), cell_number, "wheat_classic")
