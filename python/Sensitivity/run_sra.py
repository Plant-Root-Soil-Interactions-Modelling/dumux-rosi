"""
    starts a single sra simulation of soybean or maize freezing fixed parameters, 
    and passing parameters for steady state analysis
"""
import sys; sys.path.append("../modules/"); sys.path.append("../../../CPlantBox/");  sys.path.append("../../../CPlantBox/src/python_modules")
sys.path.append("../../build-cmake/cpp/python_binding/")  # DUMUX solver
import numpy as np

import scenario_setup as scenario
import evapotranspiration as evap
import sra


def run_soybean(file_name, enviro_type, sim_time, kr, kx, lmax1, lmax2, lmax3, theta1, r1, r2, a, src):

    # parameters
    soil_, table_name, min_b, max_b, cell_number, area, Kc = scenario.soybean(int(enviro_type))
    dt = 360 / (24 * 3600)  # time step [day]

    start_date = '2021-05-10 00:00:00'  # INARI csv data
    x_, y_ = evap.net_infiltration_table_beers_csvS(start_date, sim_time, evap.lai_soybean2, Kc)
    trans_soybean = evap.get_transpiration_beers_csvS(start_date, sim_time, area, evap.lai_soybean2, Kc)

    # initialize soil
    s, soil = scenario.create_soil_model(soil_, min_b, max_b, cell_number, type = 1, times = x_, net_inf = y_)
    sra_table_lookup = sra.open_sra_lookup("data/" + table_name)

    # initialize root system
    xml_name = "data/Glycine_max_Moraes2020_opt2" + "_modified" + ".xml"  # root growth model parameter file
    mods = {"lmax145":lmax1, "lmax2":lmax2, "lmax3":lmax3, "theta45":theta1, "r145":lmax1, "r2":lmax2, "r3":lmax3, "a":a, "src":src}
    r = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, xml_name, stochastic = False, mods = mods)
    # scenario.init_conductivities_const(r, kr, kx)
    # scenario.init_dynamic_simple_growth(r, kr, kr, kx, kx) # parametrisation of hydraulic conductivities
    # scenario.init_lupine_conductivities2(r, kr, kx)  # altered
    scenario.init_lupine_conductivities(r, kr, kx)

    psi_x_, psi_s_, sink_, x_, y_, psi_s2_, vol_, surf_, krs_, depth_, soil_c_, c_ = sra.simulate_dynamic(
        s, r, sra_table_lookup, sim_time, dt, trans_soybean, rs_age = 1., type_ = 1)
    # print("volume")
    # print(np.sum(vol_, axis = 0))
    # print("surface")
    # print(np.sum(surf_, axis = 0))
    # scenario.write_files(file_name, psi_x_, psi_s_, sink_, x_, y_, psi_s2_, vol_, surf_, krs_, depth_)

    np.save('results/transpiration_' + file_name, np.vstack((x_, -np.array(y_))))
    print("finished " + file_name)


def run_maize(file_name, enviro_type, sim_time, kr, kx, lmax1, lmax2, lmax3, theta1, r1, r2, a, delaySB):

    # parameters
    soil_, table_name, min_b, max_b, cell_number, area, Kc = scenario.maize(int(enviro_type))
    dt = 360 / (24 * 3600)  # time step [day]

    start_date = '2021-05-10 00:00:00'  # INARI csv data
    x_, y_ = evap.net_infiltration_table_beers_csvS(start_date, sim_time, evap.lai_maize2, Kc)
    trans_maize = evap.get_transpiration_beers_csvS(start_date, sim_time, area, evap.lai_maize2, Kc)

    # initialize soil
    s, soil = scenario.create_soil_model(soil_, min_b, max_b, cell_number, type = 2, times = x_, net_inf = y_)
    sra_table_lookup = sra.open_sra_lookup("data/" + table_name)

    xml_name = "data/Zeamays_synMRI_modified.xml"  # root growth model parameter file
    mods = {"lmax145":lmax1, "lmax2":lmax2, "lmax3":lmax3, "theta45":theta1, "r145":lmax1, "r2":lmax2, "r3":lmax3, "a":a, "delaySB":delaySB}
    r = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, xml_name, stochastic = False, mods = mods)
    # scenario.init_conductivities_const(r, kr, kx)
    # scenario.init_dynamic_simple_growth(r, kr, kr, kx, kx) # parametrisation of hydraulic conductivities
    # scenario.init_maize_conductivities2(r, kr, kx)
    # scenario.init_timing(r, kr0 = 1.e-4, kx0 = 1.e-2, dt = kr)
    scenario.init_maize_conductivities(r, kr, kx)

    psi_x_, psi_s_, sink_, x_, y_, psi_s2_, vol_, surf_, krs_, depth_, soil_c_, c_ = sra.simulate_dynamic(
        s, r, sra_table_lookup, sim_time, dt, trans_maize, rs_age = 1., type_ = 2)
    # scenario.write_files(file_name, psi_x_, psi_s_, sink_, x_, y_, psi_s2_, vol_, surf_, krs_, depth_, soil_c_, c_)
    np.save('results/transpiration_' + file_name, np.vstack((x_, -np.array(y_))))  # time [day], transpiration [cm3/day]   <- water uptake for global analysis
    np.save('results/nitrate_' + file_name, c_)  # time [day], transpiration [cm3/day]   <- water uptake for global analysis

    print("finished " + file_name)


if __name__ == "__main__":

    # print(sys.argv)
    # print(len(sys.argv[1:]))
    # file_name, enviro_type, sim_time, kr, kx, lmax0, lmax1, lmax2, theta0, r0, r1, a, src

    # file_name = sys.argv[1]
    # enviro_type = int(float(sys.argv[2]))
    # sim_time = float(sys.argv[3])
    # kr = float(sys.argv[4])
    # kx = float(sys.argv[5])
    # lmax1 = float(sys.argv[6])
    # lmax2 = float(sys.argv[7])
    # lmax3 = float(sys.argv[8])
    # theta1 = float(sys.argv[9])
    # r1 = float(sys.argv[10])
    # r2 = float(sys.argv[11])
    # a = float(sys.argv[12])
    # src = int(float(sys.argv[13]))
    # # bla

    file_name = "test_"
    enviro_type = 0
    sim_time = 1
    kr = 1.
    kx = 1.
    lmax1 = 1.
    lmax2 = 1.
    lmax3 = 1.
    theta1 = 0.
    r1 = 1.
    r2 = 1.
    a = 1.
    src = 2

    run_maize(file_name, enviro_type, sim_time, kr, kx, lmax1, lmax2, lmax3, theta1, r1, r2, a, src)

