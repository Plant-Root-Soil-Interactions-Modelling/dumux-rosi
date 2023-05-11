""" 
    Maize (water and nitrate, SRA)
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

import plantbox as pb
import visualisation.vtk_plot as vp

import scenario_setup as scenario
import evapotranspiration as evap
import sra
import functional.van_genuchten as vg

""" parameters   """
soil_, table_name, min_b, max_b, cell_number, area, Kc = scenario.maize(0)

sim_time = 95  # 47  # 95  #  [day]
dt = 360 / (24 * 3600)  # time step [day]

start_date = '2021-05-10 00:00:00'  # INARI csv data
x_, y_ = evap.net_infiltration_table_beers_csvS(start_date, sim_time, evap.lai_maize2, Kc)
trans_maize = evap.get_transpiration_beers_csvS(start_date, sim_time, area, evap.lai_maize2, Kc)

# start_date = '1995-03-15 00:00:00'
# x_, y_ = evap.net_infiltration_table_beers_pickle('data/95.pkl', start_date, sim_time, evap.lai_maize, Kc)
# trans_maize = evap.get_transpiration_beers_pickle('data/95.pkl', start_date, sim_time, area, evap.lai_maize, Kc)

""" initialize """
s, soil = scenario.create_soil_model(soil_, min_b, max_b, cell_number, type = 2, times = x_, net_inf = y_, wet = False)  # , times = x_, net_inf = y_
sra_table_lookup = sra.open_sra_lookup("data/" + table_name)
# sra_table_lookup = vg.Parameters([0.0639, 0.3698, 0.0096, 1.4646, 4.47])
# vg.create_mfp_lookup(sra_table_lookup, wilting_point = -15000, n = 15001)
# sra_table_lookup.plot_retention_curve()

xml_name = "data/Zeamays_synMRI_modified.xml"  # root growth model parameter file
r = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, xml_name)  # pass parameter file for dynamic growth
# rsml_name = "results/maize.rsml"  # created by rootsystem_maize.py
# r = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, rsml_name)
# scenario.init_conductivities_const_growth(r, kr_const = 1.8e-3, kx_const = 0.01)
# scenario.init_lupine_conductivities(r)
scenario.init_maize_conductivities(r)

""" sanity checks """
r.test()  # sanity checks

""" numerical solution """
water0 = s.getWaterVolume()  # total initial water volume in domain

psi_x_, psi_s_, sink_, x_, y_, psi_s2_, vol_, surf_, krs_, depth_, soil_c_, c_ = sra.simulate_dynamic(
    s, r, sra_table_lookup, sim_time, dt, trans_maize, rs_age = 1., type_ = 1)
# psi_x_, psi_s_, sink_, x_, y_, psi_s2_ = sra.simulate_const(s, r, sra_table_lookup, trans, sim_time, dt)

water = s.getWaterVolume()

""" output """
scenario.write_files("maize_sra0d", psi_x_, psi_s_, sink_, x_, y_, psi_s2_, vol_, surf_, krs_, depth_, soil_c_, c_)
print("\nnet water change in soil", water0 - water, "cm3")
print("fin")
