""" 
    Soybean (water only) 
"""
import sys; sys.path.append("../modules/"); sys.path.append("../../../CPlantBox/");  sys.path.append("../../../CPlantBox/src/python_modules")
sys.path.append("../../build-cmake/cpp/python_binding/")  # DUMUX solver

import scenario_setup as scenario
import evapotranspiration as evap
import sra

""" parameters   """
soil_, table_name, min_b, max_b, cell_number, area, Kc = scenario.soybean(0)

sim_time = 87.5  # 87.5  # 44  # 87.5  # [day]
dt = 360 / (24 * 3600)  # time step [day]

start_date = '2021-05-10 00:00:00'  # INARI csv data
x_, y_ = evap.net_infiltration_table_beers_csvS(start_date, sim_time, evap.lai_soybean2, Kc)
trans_soybean = evap.get_transpiration_beers_csvS(start_date, sim_time, area, evap.lai_soybean2, Kc)

# start_date = '1995-03-15 00:00:00' # pickle data from germany
# x_, y_ = evap.net_infiltration_table_beers_pickle('data/95.pkl', start_date, sim_time, evap.lai_soybean, Kc)
# trans_soybean = evap.get_transpiration_beers_pickle('data/95.pkl', start_date, sim_time, area, evap.lai_soybean, Kc)

""" initialize """
s, soil = scenario.create_soil_model(soil_, min_b, max_b, cell_number, type = 1, times = x_, net_inf = y_, wet = False)  # , times = x_, net_inf = y_
sra_table_lookup = sra.open_sra_lookup("data/" + table_name)

xml_name = "data/Glycine_max_Moraes2020_opt2_modified.xml"  # root growth model parameter file
r = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, xml_name)  # pass parameter file for dynamic growth
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

psi_x_, psi_s_, sink_, x_, y_, psi_s2_, vol_, surf_, krs_, depth_, soil_c_, c_ = sra.simulate_dynamic(
    s, r, sra_table_lookup, sim_time, dt, trans_soybean, rs_age = 1., type_ = 1)
# psi_x_, psi_s_, sink_, x_, y_, psi_s2_ = sra.simulate_const(s, r, sra_table_lookup, sim_time, dt)

water = s.getWaterVolume()

""" output """
scenario.write_files("soybean_sra0d", psi_x_, psi_s_, sink_, x_, y_, psi_s2_, vol_, surf_, krs_, depth_)
print("\nnet water change in soil", water0 - water, "cm3")
print("fin")
