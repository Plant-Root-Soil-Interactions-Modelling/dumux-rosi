""" 
    Soybean (water only) TODO
"""
import sys; sys.path.append("../modules/"); sys.path.append("../../../CPlantBox/");  sys.path.append("../../../CPlantBox/src/python_modules")
sys.path.append("../../build-cmake/cpp/python_binding/")  # DUMUX solver

from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()

import scenario_setup as scenario
from rhizo_models import RhizoMappedSegments
from xylem_flux import sinusoidal2
import cyl
import sra

""" parameters   """
min_b = [-19, -2.5, -200.]  # Domain [38 cm Reihe, 5 cm Pflanzen]
max_b = [19, 2.5, 0.]
cell_number = [1, 1, 200]  # 1 cm3

theta_r = 0.025  # sandy loam
theta_s = 0.403
alpha = 0.0383  # (cm-1) soil
n = 1.3774
k_sat = 60.  # (cm d-1)
soil_ = [theta_r, theta_s, alpha, n, k_sat]

trans = 0.6 * (38 * 5)  # cm3/day (38 * 5 = 190)

sim_time = 1.*87.5  # [day]
dt = 360 / (24 * 3600)  # time step [day]

""" initialize """
s, soil = scenario.create_soil_model(soil_, min_b, max_b, cell_number, p_top = -330, p_bot = -180, type = 1)

xml_name = "Glycine_max_Moraes2020_opt2" + "_modified" + ".xml"  # root growth model parameter file
r = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, xml_name)  # pass parameter file for dynamic growth
# rsml_name = "results/soybean.rsml"  # created by rootsystem_soybean.py
# r = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, rsml_name)

trans_f1 = lambda age, dt:-trans * sinusoidal2(age, dt) - 0.01  # Guilaumes questions
trans_f2 = lambda age, dt:-trans * sinusoidal2(age, dt) * (age / sim_time)  # growing potential transpiration

sra_table_lookup = sra.open_sra_lookup("../coupled/sra/table_jan_comp")  # make sure the soil parameters correspond to the look up table

""" sanity checks """
if rank == 0:
    r.test()  # sanity checks

""" numerical solution """
water0 = s.getWaterVolume()  # total initial water volume in domain

psi_x_, psi_s_, sink_, x_, y_, psi_s2_ = sra.simulate_dynamic(s, r, sra_table_lookup, trans, sim_time, dt, trans_f2)
# psi_x_, psi_s_, sink_, x_, y_, psi_s2_ = sra.simulate_const(s, r, sra_table_lookup, trans, sim_time, dt)

water = s.getWaterVolume()

""" output """
if rank == 0:

    scenario.write_files("soybean_sra2", psi_x_, psi_s_, sink_, x_, y_, psi_s2_)

    print("\ntotal uptake", water0 - water, "cm3")
    print("fin")
