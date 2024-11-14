"""
    starts a single sra simulation of soybean or maize freezing fixed parameters, 
    and passing parameters for steady state analysis
"""

import numpy as np
import sys

import scenario_setup as scenario

import evapotranspiration as evap
import soil_model
import hydraulic_model

import sra_new


def run_soybean(file_name, enviro_type, sim_time, kr, kx, lmax, theta1, r, a, src, kr_old = None, kx_old = None, save_all = False):
    """
        file_name                output file_name
        enviro_type              envirotype (number)
        sim_time                 simulation time [days]
        lmax                     [lmax1, lmax2, lmax3]
    """

    print("***********************")
    print("run_soybean", file_name, enviro_type, sim_time, kr, kx, lmax, theta1, r, a, src)
    print("***********************", flush = True)

    dt = 360 / (24 * 3600)  # time step [day]

    soil_, table_name, min_b, max_b, cell_number, area, Kc = scenario.soybean(int(enviro_type))  # ## TODO water table value <----
    water_table = 120.

    start_date = '2021-05-10 00:00:00'  # INARI csv data
    x_, y_ = evap.net_infiltration_table_beers_csvS(start_date, sim_time, evap.lai_soybean2, Kc)
    trans_soybean = evap.get_transpiration_beers_csvS(start_date, sim_time, area, evap.lai_soybean2, Kc)

    # initialize soil
    s = soil_model.create_richards(soil_, min_b, max_b, cell_number, times = x_, net_inf = y_, bot_bc = "potential", bot_value = 200. - water_table)
    # s = soil_model.create_richards(soil_, min_b, max_b, cell_number, times = x_, net_inf = y_, bot_bc = "noFlux", bot_value = 0.)
    print("soil model set\n", flush = True)

    # initialize root system
    print("starting hydraulic model", flush = True)
    xml_name = "data/Glycine_max_Moraes2020_opt2_modified.xml"  # root growth model parameter file

    mods = {"lmax145":lmax[0], "lmax2":lmax[1], "lmax3":lmax[2], "r145":r[0], "r2":r[1], "r3":r[1], "a":a}
    if theta1:
        mods["theta45"] = theta1
    if src:
        mods["src"] = src

    r, params = hydraulic_model.create_mapped_rootsystem(min_b, max_b, cell_number, s, xml_name, stochastic = False, mods = mods, model = "Meunier")
    print("***********************", "hydraulic model set\n", flush = True)

    if kr_old is not None:
        print("10 variable conductivity set up")  # ykr1, okr1, ykr2, okr2, kr3_, ykx1, okx1, kx2_, kx3_
        scenario.init_lupine_conductivities_sa(r, kr[0], kr_old[0], kr[1], kr_old[1], kr[2], kx[0], kx_old[0], kx[1], kx_old[1], kx[2])
    else:
        scenario.init_lupine_conductivities(r, kr, kx)
    # params.plot_conductivities(False, lateral_ind = [2, 3])  #
    # dd

    print("conductivities set\n", flush = True)

    """ sanity checks """
    r.test()  # sanity checks
    print("sanity test done\n", flush = True)

    pot_trans, psi_x_, psi_s_, sink_, x_, y_, psi_s2_, vol_, surf_, krs_, depth_, soil_c_, c_ = sra_new.simulate_dynamic(
        s, r, table_name, sim_time, dt, trans_soybean, initial_age = 1., type_ = 1)

    if save_all:
        scenario.write_results(file_name, pot_trans, psi_x_, psi_s_, sink_, x_, y_, psi_s2_, vol_, surf_, krs_, depth_)
    else:
        scenario.write_results(file_name, pot_trans, [], [], [], x_, y_, [], vol_, surf_, krs_, depth_)

    print("writing parameters", file_name)
    r.ms.writeParameters(file_name)  # corresponding xml parameter set

    # np.save('results/transpiration_' + file_name, np.vstack((x_, -np.array(y_))))
    print("finished " + file_name)


if __name__ == "__main__":

    # theta1 = None
    # src = None
    # kx = [0.1, 1.e-3, 1.e-3]
    # kx_old = [0.35, 0.015]
    # kr = [1.e-3, 4.e-3, 4.e-3]
    # kr_old = [5e-4, 0.0015]
    # run_soybean("", 0, 1, kr, kx, [1., 1., 1.], theta1, [1., 1.], 1., src, kr_old, kx_old, save_all = True)

    type = sys.argv[1]
    file_name = sys.argv[2]
    enviro_type = int(float(sys.argv[3]))
    sim_time = float(sys.argv[4])

    if type == "original":

        print("running original sa")

        kr = float(sys.argv[5])
        kx = float(sys.argv[6])
        lmax1 = float(sys.argv[7])
        lmax2 = float(sys.argv[8])
        lmax3 = float(sys.argv[9])
        theta1 = float(sys.argv[10])
        r1 = float(sys.argv[11])
        r2 = float(sys.argv[12])
        a = float(sys.argv[13])
        src = int(float(sys.argv[14]))

        run_soybean(file_name, enviro_type, sim_time, kr, kx, [lmax1, lmax2, lmax3], theta1, [r1, r2], a, src, save_all = True)  # pass only mods, kr, and kx, TODO

    elif type == "conductivities10":

        print("running conductivities 10")

        kr = np.zeros((3,))
        kr_old = np.zeros((2,))
        kx = np.zeros((3,))
        kx_old = np.zeros((2,))

        kr[0] = sys.argv[5]
        kr_old[0] = sys.argv[6]
        kr[1] = sys.argv[7]
        kr_old[1] = sys.argv[8]
        kr[2] = sys.argv[9]

        kx[0] = sys.argv[10]
        kx_old[0] = sys.argv[11]
        kx[1] = sys.argv[12]
        kx_old[1] = sys.argv[13]
        kx[2] = sys.argv[14]

        run_soybean(file_name, enviro_type, sim_time, kr, kx, [1., 1., 1.], None, [1., 1.], 1., None, kr_old, kx_old, save_all = True)

    else:

        print("Unknown run sa type")

# def run_maize(file_name, enviro_type, sim_time, kr, kx, lmax1, lmax2, lmax3, theta1, r1, r2, a, delaySB):
#
#     # parameters
#     soil_, table_name, min_b, max_b, cell_number, area, Kc = scenario.maize(int(enviro_type))
#     dt = 360 / (24 * 3600)  # time step [day]
#
#     start_date = '2021-05-10 00:00:00'  # INARI csv data
#     x_, y_ = evap.net_infiltration_table_beers_csvS(start_date, sim_time, evap.lai_maize2, Kc)
#     trans_maize = evap.get_transpiration_beers_csvS(start_date, sim_time, area, evap.lai_maize2, Kc)
#
#     # initialize soil
#     s, soil = scenario.create_soil_model(soil_, min_b, max_b, cell_number, type = 2, times = x_, net_inf = y_)
#     # sra_table_lookup = sra.open_sra_lookup("data/" + table_name) ##############################
#
#     xml_name = "data/Zeamays_synMRI_modified.xml"  # root growth model parameter file
#     mods = {"lmax145":lmax1, "lmax2":lmax2, "lmax3":lmax3, "theta45":theta1, "r145":lmax1, "r2":lmax2, "r3":lmax3, "a":a, "delaySB":delaySB}
#     r = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, xml_name, stochastic = False, mods = mods)
#     # scenario.init_conductivities_const(r, kr, kx)
#     # scenario.init_dynamic_simple_growth(r, kr, kr, kx, kx) # parametrisation of hydraulic conductivities
#     # scenario.init_maize_conductivities2(r, kr, kx)
#     # scenario.init_timing(r, kr0 = 1.e-4, kx0 = 1.e-2, dt = kr)
#     scenario.init_maize_conductivities(r, kr, kx)
#
#     psi_x_, psi_s_, sink_, x_, y_, psi_s2_, vol_, surf_, krs_, depth_, soil_c_, c_ = sra.simulate_dynamic(
#         s, r, soil, sim_time, dt, trans_maize, rs_age = 1., type_ = 2)  # soil, sra_table_lookup
#     # scenario.write_files(file_name, psi_x_, psi_s_, sink_, x_, y_, psi_s2_, vol_, surf_, krs_, depth_, soil_c_, c_)
#     np.save('results/transpiration_' + file_name, np.vstack((x_, -np.array(y_))))  # time [day], transpiration [cm3/day]   <- water uptake for global analysis
#     np.save('results/nitrate_' + file_name, c_)  # time [day], transpiration [cm3/day]   <- water uptake for global analysis
#
#     print("finished " + file_name)
