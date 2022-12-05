""" 
    Maize using rhizosphere models for a constant root system TODO part for nitrate is unstable 
"""
import sys; sys.path.append("../modules/"); sys.path.append("../../../CPlantBox/");  sys.path.append("../../../CPlantBox/src/python_modules")
sys.path.append("../../build-cmake/cpp/python_binding/")  # DUMUX solver

from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()

import scenario_setup as scenario
from rhizo_models import RhizoMappedSegments
from xylem_flux import sinusoidal2
import cyl3

"""
 CURRENTLY I am at the state to get the rhizosphere running, whitout any nitrate
"""

""" parameters   """
min_b = [-37.5, -7.5, -200.]  # Domain Mais: 60 cm Reihe, 10 cm Pflanzen
max_b = [37.5, 7.5, 0.]
cell_number = [1, 1, 200]  # 1 cm3

theta_r = 0.025  # sandy loam
theta_s = 0.403
alpha = 0.0383  # (cm-1) soil
n = 1.3774
k_sat = 60.  # (cm d-1)
soil_ = [theta_r, theta_s, alpha, n, k_sat]

trans = 0.6 * (75 * 15)  # cm3/day

sim_time = 95  #  [day]
dt = 1800 / (24 * 3600)  # time step [day] 20

""" rhizosphere model parameters """
wilting_point = -15000  # cm
nc = 10  # dof+1
logbase = 0.5  # according to Mai et al. (2019)
mode = "dumux_dirichlet"

""" initialize """
s, soil = scenario.create_soil_model(soil_, min_b, max_b, cell_number, p_top = -330, p_bot = -180, type = 2)

rsml_name = "results/maize.rsml"  # created by rootsystem_maize.py
r = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, rsml_name)

rs = RhizoMappedSegments(r, wilting_point, nc, logbase, mode)

trans_f1 = lambda age, dt:-trans * sinusoidal2(age, dt)  # Guilaumes questions - 0.01
trans_f2 = lambda age, dt:-trans * sinusoidal2(age, dt) * (age / sim_time)  # growing potential transpiration

""" sanity checks """
if rank == 0:
    r.test()  # sanity checks

""" numerical solution """
water0 = s.getWaterVolume()  # total initial water volume in domain

psi_x_, psi_s_, sink_, x_, y_, psi_s2_ = cyl3.simulate_const(s, rs, trans, sim_time, dt, trans_f2)  # trans_f

water = s.getWaterVolume()

""" output """
if rank == 0:

    scenario.write_files("maize_cyl2", psi_x_, psi_s_, sink_, x_, y_, psi_s2_)

    print("\ntotal uptake", water0 - water, "cm3")
    print("fin")
