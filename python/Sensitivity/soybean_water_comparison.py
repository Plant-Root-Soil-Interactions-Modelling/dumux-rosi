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

sim_time = 7  # a week
dt = 360 / (24 * 3600)  # time step [day]
trans_f = lambda t, dt:-HydraulicModel_Meunier.sinusoidal2(t, dt) * area * 5
soil_, table_name, min_b, max_b, cell_number, area, Kc = scenario.soybean(0)

""" 1D """
# initialize
s, soil = scenario.create_soil_model(soil_, min_b, max_b, cell_number, type = 1, times = [0, 200], net_inf = [0., 0.], bot_bc = "noFlux")  # , times = x_, net_inf = y_
xml_name = "data/Glycine_max_Moraes2020_opt2_modified.xml"  # root growth model parameter file
r = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, xml_name)  # pass parameter file for dynamic growth
scenario.init_lupine_conductivities(r)
# numerical solution """
water0 = s.getWaterVolume()  # total initial water volume in domain
psi_x_, psi_s_, sink_, times, act_trans, psi_s2_, vol_, surf_, krs_, depth_, soil_c_, c_ = sra_new.simulate_dynamic(
    s, r, "data/" + table_name, sim_time, dt, trans_f, initial_age = 10., type_ = 1)
water = s.getWaterVolume()
# output """
plot_transpiration(times, act_trans, act_trans, lambda t:-trans_f(t, dt))
scenario.write_files("soybean_comparison1D", psi_x_, psi_s_, sink_, times, act_trans, psi_s2_, vol_, surf_, krs_, depth_)
print("\nnet water change in soil 1D", water0 - water, "cm3")

# """ 2D """
# cell_number = [38, 1, 200]
# # initialize
# s, soil = scenario.create_soil_model(soil_, min_b, max_b, cell_number, type = 1, times = [0, 200], net_inf = [0., 0.], bot_bc = "noFlux")  # , times = x_, net_inf = y_
# xml_name = "data/Glycine_max_Moraes2020_opt2_modified.xml"  # root growth model parameter file
# r = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, xml_name)  # pass parameter file for dynamic growth
# scenario.init_lupine_conductivities(r)
# # numerical solution """
# water0 = s.getWaterVolume()  # total initial water volume in domain
# psi_x_, psi_s_, sink_, times, act_trans, psi_s2_, vol_, surf_, krs_, depth_, soil_c_, c_ = sra_new.simulate_dynamic(
#     s, r, "data/" + table_name, sim_time, dt, trans_f, initial_age = 10., type_ = 1)
# water = s.getWaterVolume()
# # output """
# plot_transpiration(times, act_trans, act_trans, lambda t:-trans_f(t, dt))
# scenario.write_files("soybean_comparison1D", psi_x_, psi_s_, sink_, times, act_trans, psi_s2_, vol_, surf_, krs_, depth_)
# print("\nnet water change in soil 1D", water0 - water, "cm3")

print("fin")
