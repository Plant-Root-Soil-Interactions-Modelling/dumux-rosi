"""
    starts a single sra simulation of soybean or maize freezing fixed parameters, 
    and passing parameters for steady state analysis
"""
import sys; sys.path.append("../modules/"); sys.path.append("../../../CPlantBox/");  sys.path.append("../../../CPlantBox/src/python_modules")
sys.path.append("../../build-cmake/cpp/python_binding/")  # DUMUX solver

import scenario_setup as scenario
import evapotranspiration as evap
import sra


def run_soybean(file_name, enviro_type, sim_time, kr, kx, lmax1, lmax2, lmax3, theta1, r1, r2, a, src):

    # parameters
    soil_, table_name, p_top, min_b, max_b, cell_number, area, Kc = scenario.maize(0)

    sim_time = 1  # 87.5  # [day]
    dt = 360 / (24 * 3600)  # time step [day]

    start_date = '1995-03-15 00:00:00'
    x_, y_ = evap.net_infiltration_table_beers('data/95.pkl', start_date, sim_time, evap.lai_soybean, Kc)
    trans_soybean = evap.get_transpiration_beers('data/95.pkl', start_date, sim_time, area, evap.lai_soybean, Kc)

    # initialize soil
    s, soil = scenario.create_soil_model(soil_, min_b, max_b, cell_number, p_top = p_top, p_bot = (p_top + 200), type = 1, times = x_, net_inf = y_)
    sra_table_lookup = sra.open_sra_lookup("data/" + table_name)

    # initialize root system
    xml_name = "Glycine_max_Moraes2020_opt2" + "_modified" + ".xml"  # root growth model parameter file
    mods = {"lmax145":lmax1, "lmax2":lmax2, "lmax3":lmax3, "theta45":theta1, "r145":r1, "r2":r2, "a":a, "src":src}
    r = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, xml_name, stochastic = False, mods = mods)
    # scenario.init_lupine_conductivities(r)
    scenario.init_dynamic_simple_growth(r, kr, kr, kx, kx)
    # scenario.init_conductivities_const(r, kr, kx)

    # numerical solution
    psi_x_, psi_s_, sink_, x_, y_, psi_s2_, vol_, surf_, krs_, depth_, soil_c_, c_ = sra.simulate_dynamic(s, r, sra_table_lookup, sim_time, dt, trans_soybean)
    scenario.write_files(file_name, psi_x_, psi_s_, sink_, x_, y_, psi_s2_, vol_, surf_, krs_, depth_)
    print("finished " + file_name)


if __name__ == "__main__":

    # print(sys.argv)
    # print(len(sys.argv[1:]))
    # file_name, enviro_type, sim_time, kr, kx, lmax0, lmax1, lmax2, theta0, r0, r1, a, src

    file_name = sys.argv[1]
    enviro_type = int(float(sys.argv[2]))
    sim_time = float(sys.argv[3])
    kr = float(sys.argv[4])
    kx = float(sys.argv[5])
    print("kr", kr, "kx", kx)
    lmax1 = float(sys.argv[6])
    lmax2 = float(sys.argv[7])
    lmax3 = float(sys.argv[8])
    theta1 = float(sys.argv[9])
    r1 = float(sys.argv[10])
    r2 = float(sys.argv[11])
    a = float(sys.argv[12])
    src = int(float(sys.argv[13]))

    run_soybean(file_name, enviro_type, sim_time, kr, kx, lmax1, lmax2, lmax3, theta1, r1, r2, a, src)

