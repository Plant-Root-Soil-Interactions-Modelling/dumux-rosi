""" 
Single root scenario - soil depletion due to sinusoidal transpiration over 21 days

using steady rate approach and fix point iteration (sra)
"""
import sys; sys.path.append("../modules/"); sys.path.append("../../../CPlantBox/");  sys.path.append("../../../CPlantBox/src/python_modules")

import scenario_setup as scenario
import sra

""" parameters   """
min_b = [-6, -1.5, -150.]  # domain 12cm x 3cm x 150cm
max_b = [6, 1.5, 0.]
cell_number = [1, 1, 150]  # 12x3x1 cm3 - 1d soil

theta_r = 0.025
theta_s = 0.403
alpha = 0.0383  # (cm-1) soil
n = 1.3774
k_sat = 60.  # (cm d-1)
soil_ = [theta_r, theta_s, alpha, n, k_sat]

trans = 0.6 * (12 * 3)  # cm3/day

sim_time = 21  #  [day]
dt = 60 / (24 * 3600)  # time step [day]

""" initialize """
s, soil = scenario.create_soil_model(soil_, min_b, max_b, cell_number, p_top = -330, p_bot = -180)
r = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, "results/wheat_nobasals.rsml")  # created by rootsystem.py
sra_table_lookup = sra.open_sra_lookup("../coupled/sra/table_jan_comp")  # make sure the soil parameters correspond to the look up table
# sra_table_lookup = soil  # without using the lookup table.

""" sanity checks """
r.test()  # we might add more
# print("Krs", r.get_krs(0.))

""" numerical solution """
water0 = s.getWaterVolume()  # total initial water volume in domain

psi_x_, psi_s_, sink_, x_, y_, psi_s2_ = sra.simulate_const(s, r, sra_table_lookup, trans, sim_time, dt)

scenario.write_files("rootsystem_nobasals1d_sra", psi_x_, psi_s_, sink_, x_, y_, psi_s2_)  # 1d soil

print("\ntotal uptake", water0 - s.getWaterVolume(), "cm3")
print("fin")

