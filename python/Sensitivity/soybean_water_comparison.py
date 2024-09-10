""" 
    Simulates water movement for a rather simple scenario and compares errors 2D vs 1D 
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

import matplotlib
import matplotlib.pyplot as plt

import functional.van_genuchten as vg
from functional.PlantHydraulicModel import HydraulicModel_Meunier

from rhizo_models import plot_transpiration

import scenario_setup as scenario
import evapotranspiration as evap
import sra_new

sim_time = 1  # a week
dt = 360 / (24 * 3600)  # time step [day]
trans_f = lambda t, dt:-HydraulicModel_Meunier.sinusoidal2(t, dt) * area * 5
soil_, table_name, min_b, max_b, cell_number, area, Kc = scenario.soybean(0)

initial_age = 10

""" 1D """
# initialize
s, soil = scenario.create_soil_model(soil_, min_b, max_b, cell_number, type = 1, times = [0, 200], net_inf = [0., 0.], bot_bc = "noFlux")  # , times = x_, net_inf = y_
xml_name = "data/Glycine_max_Moraes2020_opt2_modified.xml"  # root growth model parameter file
r = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, xml_name)  # pass parameter file for dynamic growth
scenario.init_lupine_conductivities(r)
r.ms.simulate(initial_age)
# numerical solution """
water0 = s.getWaterVolume()  # total initial water volume in domain
psi_x_, psi_s_, sink_, times, act_trans, psi_s2_, vol_, surf_, krs_, depth_, soil_c_, c_ = sra_new.simulate_dynamic(
    s, r, "data/" + table_name, sim_time, dt, trans_f, initial_age = initial_age, type_ = 1)
water = s.getWaterVolume()
# output """
scenario.write_files("soybean_comparison1D", psi_x_, psi_s_, sink_, times, act_trans, psi_s2_, vol_, surf_, krs_, depth_)
print("\nnet water change in soil 1D", water0 - water, "cm3")
net_change1d = water0 - water

# """ 3D """
# cell_number = [38, 3, 200]
# # initialize
# s, soil = scenario.create_soil_model(soil_, min_b, max_b, cell_number, type = 1, times = [0, 200], net_inf = [0., 0.], bot_bc = "noFlux")  # , times = x_, net_inf = y_
# xml_name = "data/Glycine_max_Moraes2020_opt2_modified.xml"  # root growth model parameter file
# r = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, xml_name)  # pass parameter file for dynamic growth
# scenario.init_lupine_conductivities(r)
# r.ms.simulate(initial_age)
# picker = lambda x, y, z: s.pick([x, y, z])  #  function that return the index of a given position in the soil grid (should work for any grid - needs testing)
# r.ms.setSoilGrid(picker)  # maps segments, maps root segements and soil grid indices to each other in both directions
#
# # numerical solution """
# water0 = s.getWaterVolume()  # total initial water volume in domain
# psi_x_, psi_s_, sink_, times2d, act_trans2d, psi_s2_, vol_, surf_, krs_, depth_, soil_c_, c_ = sra_new.simulate_dynamic(
#     s, r, "data/" + table_name, sim_time, dt, trans_f, initial_age = initial_age, type_ = 1)
# water = s.getWaterVolume()
# # output """
# scenario.write_files("soybean_comparison2D", psi_x_, psi_s_, sink_, times2d, act_trans2d, psi_s2_, vol_, surf_, krs_, depth_)
# print("\nnet water change in soil 2D", water0 - water, "cm3")
# net_change2d = water0 - water

print("fin")

print("net_change 1d", net_change1d)
print("net_change 2d", net_change2d)

plot_transpiration(times, act_trans, act_trans, lambda t:-trans_f(t, dt))
plot_transpiration(times2d, act_trans2d, act_trans2d, lambda t:-trans_f(t, dt))
