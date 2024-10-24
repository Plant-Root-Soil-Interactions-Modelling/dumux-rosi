""" 
    Simulates water movement for a single fixed soybean scenario
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

import functional.van_genuchten as vg

import scenario_setup as scenario

import evapotranspiration as evap
import soil_model
import hydraulic_model

import sra
import sra_new

""" parameters   """
sim_time = 87.5  # [day]
dt = 360 / (24 * 3600)  # time step [day]

soil_, table_name, min_b, max_b, cell_number, area, Kc = scenario.soybean(0)

start_date = '2021-05-10 00:00:00'  # INARI csv data
x_, y_ = evap.net_infiltration_table_beers_csvS(start_date, sim_time, evap.lai_soybean2, Kc)
trans_soybean = evap.get_transpiration_beers_csvS(start_date, sim_time, area, evap.lai_soybean2, Kc)

# start_date = '1995-03-15 00:00:00' # pickle data from germany
# x_, y_ = evap.net_infiltration_table_beers_pickle('data/95.pkl', start_date, sim_time, evap.lai_soybean, Kc)
# trans_soybean = evap.get_transpiration_beers_pickle('data/95.pkl', start_date, sim_time, area, evap.lai_soybean, Kc)

""" initialize """
# s, soil = scenario.create_soil_model(soil_, min_b, max_b, cell_number, type = 1, times = x_, net_inf = y_)  # , times = x_, net_inf = y_
s = soil_model.create_richards(soil_, min_b, max_b, cell_number, times = x_, net_inf = y_, bot_bc = "noFlux", bot_value = 0.)

xml_name = "data/Glycine_max_Moraes2020_opt2_modified.xml"  # root growth model parameter file
r = hydraulic_model.create_mapped_rootsystem(min_b, max_b, cell_number, s, xml_name, stochastic = False, mods = None, model = "Doussan")  # pass parameter file for dynamic growth

# picker = lambda x, y, z: s.pick([0., 0., z])  #  function that return the index of a given position in the soil grid (should work for any grid - needs testing)
# r.ms.setSoilGrid(picker)

# scenario.init_conductivities_const_growth(r)
scenario.init_lupine_conductivities(r)
# scenario.init_dynamic_simple_growth(r, 1.e-3, 4.e-3, 5.e-2, 2.e-3)
# r.plot_conductivities(monocot = True, plot_now = True, axes_ind = [1, 4, 5], lateral_ind = [2, 3])
# rsml_name = "results/soybean.rsml"  # created by rootsystem_soybean.py
# r = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, rsml_name)

""" sanity checks """
r.test()  # sanity checks

""" numerical solution """
water0 = s.getWaterVolume()  # total initial water volume in domain

psi_x_, psi_s_, sink_, x_, y_, psi_s2_, vol_, surf_, krs_, depth_, soil_c_, c_ = sra_new.simulate_dynamic(
    s, r, table_name, sim_time, dt, trans_soybean, 1., type_ = 1)  # , model = "Doussan"

water = s.getWaterVolume()

""" output """
scenario.write_files("soybean_noFlux", psi_x_, psi_s_, sink_, x_, y_, psi_s2_, vol_, surf_, krs_, depth_)
print("\nnet water change in soil", water0 - water, "cm3")
print("fin")
