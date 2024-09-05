""" 
    Simulates water movement for a single fixed soybean scenario
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

import functional.van_genuchten as vg

import scenario_setup as scenario
import evapotranspiration as evap
import sra_new

""" parameters   """
soil_, table_name, min_b, max_b, cell_number, area, Kc = scenario.soybean(0)

sim_time = 87.5  # [day]
dt = 360 / (24 * 3600)  # time step [day]

start_date = '2021-05-10 00:00:00'  # INARI csv data
x_, y_ = evap.net_infiltration_table_beers_csvS(start_date, sim_time, evap.lai_soybean2, Kc)
trans_soybean = evap.get_transpiration_beers_csvS(start_date, sim_time, area, evap.lai_soybean2, Kc)

""" initialize """
s, soil = scenario.create_soil_model(soil_, min_b, max_b, cell_number, type = 1, times = x_, net_inf = y_, wet = False)  # , times = x_, net_inf = y_

xml_name = "data/Glycine_max_Moraes2020_opt2_modified.xml"  # root growth model parameter file
r = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, xml_name)  # pass parameter file for dynamic growth

scenario.init_lupine_conductivities(r)

""" sanity checks """
r.test()  # sanity checks

""" numerical solution """
water0 = s.getWaterVolume()  # total initial water volume in domain

psi_x_, psi_s_, sink_, x_, y_, psi_s2_, vol_, surf_, krs_, depth_, soil_c_, c_ = sra_new.simulate_dynamic(
    s, r, "data/" + table_name, sim_time, dt, trans_soybean, initial_age = 1., type_ = 1)

water = s.getWaterVolume()

""" output """
scenario.write_files("soybean", psi_x_, psi_s_, sink_, x_, y_, psi_s2_, vol_, surf_, krs_, depth_)
print("\nnet water change in soil", water0 - water, "cm3")
print("fin")
