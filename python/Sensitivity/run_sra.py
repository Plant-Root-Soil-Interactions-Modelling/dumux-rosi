"""
    Dynamic

    Starts a single dynamic simulation (dynamic water movement with nonlinear conductivities in the perirhizal zone), by calling run_soybean() 
        
    __main__ takes command line arguments, e.g. produced by run_SA.py (to run multiple simulations on the cluster)
    interpretation of the arguments is based on type = sys.argv[1], 
    and passed to run_soybean(exp_name, enviro_type, sim_time, mods, save_all = False), 
    where the dictionary mods holds most parameters
    
    see also sra_new.simulate_dynamic(..) for the simulation loop
        
    see also run_cplantbox.py (for the macroscopic simulation)
    
    Daniel Leitner, 2025  
"""

import sys
import os
import scipy
import copy
import zipfile
import json
import numpy as np
import matplotlib.pyplot as plt

import scenario_setup as scenario
import evapotranspiration as evap
import soil_model
import hydraulic_model
import sra_new


class NpEncoder(json.JSONEncoder):

    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.bool):
            return bool(obj)
        return super(NpEncoder, self).default(obj)


def run_soybean_exp(exp_name, enviro_type, sim_time, save_all = False):
    """
    reproduces a simulation with parameters given in file "soybean_all14.zip"
    """
    folder_path = "results_cplantbox/"
    file_path = os.path.join(folder_path, exp_name + "_mods.json")
    with open(file_path, 'r', encoding = 'utf-8') as file:
        params = json.load(file)

    params["mecha_path"] = "mecha_results"  # mecha results must correspond to conductivity indices given in soybean_all14.zip and results_cplantbox/
    assert exp_name == params["exp_name"], "run_sra() type == 'file': something is wrong with exp_name"
    params.pop("exp_name")
    params.pop("enviro_type")
    params.pop("sim_time")
    run_soybean(exp_name, enviro_type, sim_time, params, save_all = True)


def run_soybean(exp_name, enviro_type, sim_time, mods, save_all = False):
    """
        exp_name                experiment name (name for output files)
        enviro_type             envirotype (number)
        sim_time                simulation time [days]
        mods                    parameters that are adjusted from the base parametrization (default: data/Glycine_max_Moraes2020_opt2_modified.xml)
        save_all                ---
    """

    if not "conductivity_mode" in mods:
            print("run_sra.run_soybean() conductivity_mode and output_times must be defined in mods dictionary")
            raise()

    if not "output_times" in mods:
        mods["output_times"] = [sim_time / 2]  # add mid for now

    sim_params = {"exp_name": exp_name, "enviro_type": enviro_type, "sim_time":sim_time}

    print("***********************")
    print("run_sra.run_soybean", exp_name, enviro_type, sim_time)
    print(mods)
    print("***********************", flush = True)
    mods_copy = copy.deepcopy(mods)
    mods_copy.update(sim_params)  # all input arguments, to make reproducable

    if "filename" in mods:
        xml_name = mods["filename"]
        mods.pop("filename")
    else:
        xml_name = "data/Glycine_max_Moraes2020_opt2_modified.xml"  # root growth model parameter file

    if "initial_age" in mods:
        initial_age = mods["initial_age"]
    else:
        initial_age = 1.  # day

    if "initial_totalpotential" in mods:
        initial_totalpotential = mods["initial_totalpotential"]
        mods.pop("initial_totalpotential")
    else:
        initial_totalpotential = -100.  # cm

    if "water_table" in mods:
        water_table = mods["water_table"]
        mods.pop("water_table")
    else:
        water_table = 120.  # cm at -200cm

    dt = 360 / (24 * 3600)  # time step [day] ####################################### currently hard coded

    soil_, table_name, min_b, max_b, cell_number, area, Kc = scenario.soybean(int(enviro_type))

    if "domain_size" in mods:
        domain_size = np.array(mods["domain_size"])
        min_b = -np.array([domain_size[0] / 2., domain_size[1] / 2., domain_size[2]])
        max_b = np.array([domain_size[0] / 2., domain_size[1] / 2., 0.])
        area = domain_size[0] * domain_size[1]
        # print("Changed domain size to", min_b, max_b, area)
        mods.pop("domain_size")

    start_date = '2021-05-10 00:00:00'  # INARI csv data
    x_, y_ = evap.net_infiltration_table_beers_csvS(start_date, sim_time, evap.lai_soybean2, Kc, initial_age)
    trans_soybean = evap.get_transpiration_beers_csvS(start_date, sim_time, area, evap.lai_soybean2, Kc, initial_age)  # lambda t, d
    if "scale_trans" in mods:
        scale_trans = mods["scale_trans"]
        trans_ = lambda t, dt: trans_soybean(t, dt) * scale_trans
        mods.pop("scale_trans")
    else:
        trans_ = trans_soybean

    if "bot_bc" in mods:
        bot_bc = mods["bot_bc"]
        mods.pop("bot_bc")
    else:
        bot_bc = "potential"

    s = soil_model.create_richards(soil_, min_b, max_b, cell_number, times = x_, net_inf = y_,
                                   bot_bc = bot_bc, bot_value = 200. - water_table, initial_totalpotential = initial_totalpotential)
    water0 = s.getWaterVolume()
    # print("soil model set\n", flush = True)

    # initialize root system
    # print("starting hydraulic model", flush = True)
    cdata = scenario.prepare_conductivities(mods)
    r, params = hydraulic_model.create_mapped_rootsystem(min_b, max_b, cell_number, s, xml_name, stochastic = False, mods = mods, model = "Doussan")  # Meunier
    scenario.set_conductivities(params, mods, cdata)
    # print("***********************", "hydraulic model set\n", flush = True)

    output_times = mods["output_times"]  # TODO what to do with it in a dynamic context
    mods.pop("output_times")

    if mods:  # something unused in mods
        print("\n********************************************************************************")
        print("scenario_setup.create_mapped_rootsystem() WARNING unused parameters:")
        for k, v in mods.items():
            print("key:", k)
        print("********************************************************************************\n")
        # raise

    # if kr_old is not None:
    #     print("10 variable conductivity set up")  # ykr1, okr1, ykr2, okr2, kr3_, ykx1, okx1, kx2_, kx3_
    #     # print(kr)
    #     # print(kr_old)
    #     # print(kr[0], kr[1])
    #     scenario.init_lupine_conductivities_sa(r.params, kr[0], kr_old[0], kr[1], kr_old[1], kr[2], kx[0], kx_old[0], kx[1], kx_old[1], kx[2])
    # else:
    #     scenario.init_lupine_conductivities(r.params, kr, kx)

    # params.plot_conductivities(False, lateral_ind = [2, 3])  #
    # dd
    # print("conductivities set\n", flush = True)

    """ sanity checks """
    r.test()  # sanity checks
    # print("sanity test done\n", flush = True)
    try:
        pot_trans, psi_x_, psi_s_, sink_, x_, y_, psi_s2_, vol_, surf_, krs_, depth_, collar_pot_ = sra_new.simulate_dynamic(
            s, r, table_name, sim_time, dt, trans_, initial_age = initial_age)
    except:
        print("***********************")
        print("EXCEPTION run_soybean", exp_name, enviro_type, sim_time)
        print(mods)
        print("***********************", flush = True)
        raise

    print("writing", exp_name + "_" + str(enviro_type))
    if save_all:
        scenario.write_results(exp_name + "_" + str(enviro_type), pot_trans, psi_x_, psi_s_, sink_, x_, y_, psi_s2_, vol_, surf_, krs_, depth_, collar_pot_)
    else:
        pass
        # scenario.write_results(exp_name, pot_trans, [], [], [], x_, y_, [], vol_, surf_, krs_, depth_)
    with open("results/" + exp_name + "_" + str(enviro_type) + "_mods.json", "w+") as f:
        json.dump(mods_copy, f, cls = NpEncoder)
    r.ms.writeParameters("results/" + exp_name + ".xml")  # corresponding xml parameter set

    print("finished " + exp_name)

    cu = scipy.integrate.simpson(y_, x = x_)  # cumulative uptake
    mean_pot = np.mean(collar_pot_)  # mean collar potential

    water_end = s.getWaterVolume()

    print("cumulative uptake is", cu, np.sum(np.multiply(y_[1:], np.diff(x_))))
    print("net water change is", water_end - water0)

    return cu


if __name__ == "__main__":

    # enviro_type = 0
    # sim_time = 1
    # theta1 = None
    # src = None
    # # kx = [0.1, 1.e-3, 1.e-3]  # cm3/day
    # # kx_old = [0.35, 0.015]
    # # kr = [1.e-3, 4.e-3, 4.e-3]  # 1/day
    # # kr_old = [5e-4, 0.0015]
    # mods = {
    #     "output_times": [40],
    #     "conductivity_mode": "scale",
    #     "scale_kr":1.,
    #     "scale_kx":1.,
    #     }
    # run_soybean("test", enviro_type, sim_time, mods, save_all = True)

    type = sys.argv[1]
    exp_name = sys.argv[2]
    enviro_type = int(float(sys.argv[3]))
    sim_time = float(sys.argv[4])

    if type == "original":  # relevant parameters
        print("running original sa", flush = True)
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
        mods = {
            "conductivity_mode": "scale",
            "scale_kr": kr,
            "scale_kx": kx,
            "lmax145": lmax1, "lmax2": lmax2, "lmax3": lmax3,
            "r145": r1, "r2": r2, "r3": r2,
            "a":a,
            "theta45": theta1,
            "src": src
            }
        run_soybean(exp_name, enviro_type, sim_time, mods, save_all = True)
    elif type == "original_new2":
        print("running new2", flush = True)  # ["r", "r145", "r2", "r3", "ln", "ln145", "ln2",  "a145", "a2", "a3"]
        r = float(sys.argv[5])
        r145 = float(sys.argv[6])
        r2 = float(sys.argv[7])
        r3 = float(sys.argv[8])
        ln = float(sys.argv[9])
        ln145 = float(sys.argv[10])
        ln2 = float(sys.argv[11])
        a145 = float(sys.argv[12])
        a2 = float(sys.argv[13])
        a3 = int(float(sys.argv[14]))
        mods = {
            "conductivity_mode": "scale",
            "r": r, "r145": r145, "r2": r2, "r3": r3,
            "ln": ln, "ln145": ln145, "ln2": ln2,
            "a145":a145, "a2": a2, "a3": a3
            }
        run_soybean(exp_name, enviro_type, sim_time, mods, save_all = True)
    elif type == "conductivities10":  # makes little sense to treat kr, kx independently
        print("running conductivities10", flush = True)
        mods = {
            "conductivity_mode": "age_dependent",
            "kr_young1": float(sys.argv[5]),
            "kr_old1": float(sys.argv[6]),
            "kr_young2": float(sys.argv[7]),
            "kr_old2": float(sys.argv[8]),
            "kr3": float(sys.argv[9]),
            "kx_young1": float(sys.argv[10]),
            "kx_old1": float(sys.argv[11]),
            "kx_young2": float(sys.argv[12]),
            "kx_old2": float(sys.argv[13]),
            "kx3": float(sys.argv[14])
        }
        run_soybean(exp_name, enviro_type, sim_time, mods, save_all = True)

    elif type == "singleroot_conductivities10":  # a single root example
        print("singleroot_conductivities10", flush = True)
        mods = {
            "conductivity_mode": "age_dependent",
            "kr_young1": float(sys.argv[5]),
            "kr_old1": float(sys.argv[6]),
            "kr_young2": float(sys.argv[5]),
            "kr_old2": float(sys.argv[6]),
            "kr3": float(sys.argv[5]),
            "kx_young1": float(sys.argv[7]),
            "kx_old1": float(sys.argv[8]),
            "kx_young2": float(sys.argv[7]),
            "kx_old2": float(sys.argv[8]),
            "kx3": float(sys.argv[7]),
            "filename": "data/Glycine_max_Moraes2020_singleroot.xml",
            "initial_age": 1.,
            "initial_totalpotential":-1000,
            "domain_size": [1., 1., 200.],
            "bot_bc": "noFlux",
        }
        run_soybean(exp_name, enviro_type, sim_time, mods, save_all = True)

    elif type == "tropisms":
        print("running tropisms", flush = True)
        n145, n2, n3 = float(sys.argv[5]), float(sys.argv[6]), float(sys.argv[7])
        s145, s2, s3 = float(sys.argv[8]), float(sys.argv[9]), float(sys.argv[10])
        mods = {
                "conductivity_mode": "scale",
                "scale_kr":1.,
                "scale_kx":1.,
                "tropismN145":n145, "tropismN2":n2, "tropismN3":n3,
                "tropismS145":s145, "tropismS2":s2, "tropismS3":s3
                }
        run_soybean(exp_name, enviro_type, sim_time, mods, save_all = True)

    elif type == "radii":
        print("running radii", flush = True)
        a145, a2, a3 = float(sys.argv[5]), float(sys.argv[6]), float(sys.argv[7])
        hairsZone145, hairsZone2, hairsZone3 = float(sys.argv[8]), float(sys.argv[9]), float(sys.argv[10])
        hairsLength145, hairsLength2, hairsLength3 = float(sys.argv[11]), float(sys.argv[12]), float(sys.argv[13])
        mods = {
            "conductivity_mode": "scale",
            "scale_kr":1.,
            "scale_kx":1.,
            "a145":a145, "a":a2, "a3":a3,
            "hairsZone145":hairsZone145, "hairsZone2":hairsZone145, "hairsZone3":hairsZone145,
            "hairsLength145":hairsLength145, "hairsLength2":hairsLength2, "hairsLength3":hairsLength3
            }

        run_soybean(exp_name, enviro_type, sim_time, mods, save_all = True)

    elif type == "file":

        with zipfile.ZipFile("soybean_all14.zip", "r") as zipf:  # hard coded if zip files are used
            with zipf.open("soybean_all14.json", "r") as json_file:
                all = json.load(json_file)  # Deserialize JSON data

        # folder_path = "results_cplantbox/"
        # file_path = os.path.join(folder_path, exp_name + "_mods.json")
        # with open(file_path, 'r', encoding = 'utf-8') as file:
        #     params = json.load(file)

        params = all[exp_name]
        params["mecha_path"] = "mecha_results"  # mecha results must correspond to conductivity indices given in soybean_all14.zip

        assert exp_name == params["exp_name"], "run_sra() type == 'file': something is wrong with exp_name"
        params.pop("exp_name")
        params.pop("enviro_type")
        params.pop("sim_time")

        run_soybean(exp_name, enviro_type, sim_time, params, save_all = True)

    else:

        print("Unknown run sa type")

# def run_maize(exp_name, enviro_type, sim_time, kr, kx, lmax1, lmax2, lmax3, theta1, r1, r2, a, delaySB):
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
#     # scenario.write_files(exp_name, psi_x_, psi_s_, sink_, x_, y_, psi_s2_, vol_, surf_, krs_, depth_, soil_c_, c_)
#     np.save('results/transpiration_' + exp_name, np.vstack((x_, -np.array(y_))))  # time [day], transpiration [cm3/day]   <- water uptake for global analysis
#     np.save('results/nitrate_' + exp_name, c_)  # time [day], transpiration [cm3/day]   <- water uptake for global analysis
#
#     print("finished " + exp_name)
