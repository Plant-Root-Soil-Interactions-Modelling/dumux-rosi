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


def run_soybean_exp(exp_name, enviro_type, sim_time, result_path, save_all = True):
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
    run_soybean(exp_name, enviro_type, sim_time, params, result_path, save_all)


def run_soybean(exp_name, enviro_type, sim_time, mods, result_path, save_all = True):
    """
        exp_name                experiment name (name for output files)
        enviro_type             envirotype (number)
        sim_time                simulation time [days]
        mods                    parameters that are adjusted from the base parametrization (default: data/Glycine_max_Moraes2020_opt2_modified.xml)
        result_path             path folder for exporting results
        save_all                additionally saves all root geometry and hydraulic data at start and end of simulation, and at time defined in mods["output_times"]
    """

    if not "conductivity_mode" in mods:
            print("run_sra.run_soybean() conductivity_mode must be defined in mods dictionary")
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
    if mods["conductivity_mode"] == "from_mecha":  # untested
        anatomy_ind = scenario.sort_indices(mods)
        anatomy = [scenario_setup.get_anatomy(ind[0]) for i in anatomy_ind]
        radii = [ mods["a145_a"], mods["a2_a"], mods["a3_a"] ]
        assert data[0, 2, 2] == radii[0], "radii failed for typ 145"  # untested (can be removed later...)
        assert data[1, 2, 2] == radii[1], "radii failed for typ 2"
        assert data[2, 2, 2] == radii[2], "radii failed for typ 3"
        cc_data = { "radii": radii, "anatomy": anatomy}
    else:
        cc_data = None
    r, hparams = hydraulic_model.create_mapped_rootsystem(min_b, max_b, cell_number, s, xml_name, stochastic = False, mods = mods, model = "Doussan")  # Meunier
    scenario.set_conductivities(hparams, mods, cdata)
    # print("***********************", "hydraulic model set\n", flush = True)

    output_times = mods["output_times"]
    output_times.sort
    mods.pop("output_times")

    if mods:  # something unused in mods
        print("\n********************************************************************************")
        print("scenario_setup.create_mapped_rootsystem() WARNING unused parameters:")
        for k, v in mods.items():
            print("key:", k)
        print("********************************************************************************\n")

    """ sanity checks """
    # r.test()  # sanity checks
    # print("sanity test done\n", flush = True)

    try:
        r1, r2, r3 = sra_new.simulate_dynamic(s, r, table_name, sim_time, dt, trans_, output_times,
                                              initial_age = initial_age, cc_data = cc_data)
    except:
        print("***********************")
        print("EXCEPTION run_soybean", exp_name, enviro_type, sim_time)
        print("***********************", flush = True)
        raise

    water_end = s.getWaterVolume()

    # write simulation results
    print("writing", exp_name + "_" + str(enviro_type))

    if not os.path.exists(result_path):
        os.mkdir(result_path)

    if save_all:
        scenario.write_results(result_path + exp_name, r1, r2, r3)
    else:
        scenario.write_results(result_path + exp_name, r1, r2, None)

    with open(result_path + exp_name + "_mods.json", "w+") as f:  # write corresponding input parameters
        json.dump(mods_copy, f, cls = NpEncoder)

    r.ms.writeParameters(result_path + exp_name + ".xml")  # corresponding xml parameter set

    # possible objective for optimization
    times = r1[0]
    act_trans = r1[2]
    collar_pot = r1[3]
    cu = scipy.integrate.simpson(act_trans, x = times)  # cumulative uptake
    mean_pot = np.mean(collar_pot)  # mean collar potential

    print("finished " + exp_name)
    print("cumulative uptake is", cu, np.sum(np.multiply(act_trans[1:], np.diff(times))))
    print("net water change is", water_end - water0)

    return cu


def run(argv):

    type = argv[1]
    exp_name = argv[2]  # codes for bot BC "_free" and "_200"
    enviro_type = int(float(argv[3]))
    sim_time = float(argv[4])

    if type == "original":  # relevant parameters
        print("running original sa", flush = True)
        kr, kx = float(argv[5]), float(argv[6])
        lmax1, lmax2, lmax3 = float(argv[7]), float(argv[8]), float(argv[9])
        theta1 = float(argv[10])
        r1, r2 = float(argv[11]), float(argv[12])  # 1., 1. -> moved to running original2
        a = float(argv[13])
        src = int(float(argv[14]))
        mods = {
            "conductivity_mode": "scale",
            "scale_kr": kr, "scale_kx": kx,
            "lmax145": lmax1, "lmax2": lmax2, "lmax3": lmax3,
            "theta45_a": theta1,
            "r145": r1, "r2": r2, "r3": r2,
            "a":a,
            "src_a": src
            }
        if "_free" in exp_name:
            mods["bot_bc"] = "freeDrainage"  # otherwise potential
        if "_200" in exp_name:
            mods["water_table"] = 200  # otherwise 120
        run_soybean(exp_name, enviro_type, sim_time, mods, exp_name[:-4] + "/", save_all = True)  # result filename = exp_name + _envirotype

    elif type == "original2":
        print("running original2", flush = True)  # ["r", "r145", "r2", "r3", "ln", "ln145", "ln2",  "a145", "a2", "a3"]
        r, r145, r2, r3 = float(argv[5]), float(argv[6]), float(argv[7]), float(argv[8])
        ln, ln145, ln2 = float(argv[9]), float(argv[10]), float(argv[11])
        a145, a2, a3 = float(argv[12]), float(argv[13]), float(argv[14])
        mods = {
            "conductivity_mode": "scale",
            "scale_kr":1., "scale_kx":1.,
            "r": r, "r145": r145, "r2": r2, "r3": r3,
            "ln": ln, "ln145": ln145, "ln2": ln2,
            "a145":a145, "a2": a2, "a3": a3
            }
        if "_free" in exp_name:
            mods["bot_bc"] = "freeDrainage"  # otherwise potential
        if "_200" in exp_name:
            mods["water_table"] = 200  # otherwise 120
        run_soybean(exp_name, enviro_type, sim_time, mods, exp_name[:-4] + "/", save_all = True)

    elif type == "tropisms":
        print("running tropisms", flush = True)
        n145, n2, n3 = float(argv[5]), float(argv[6]), float(argv[7])
        s145, s2, s3 = float(argv[8]), float(argv[9]), float(argv[10])
        theta2, theta3 = float(argv[11]), float(argv[12])
        mods = {
                "conductivity_mode": "scale",
                "scale_kr":1., "scale_kx":1.,
                "tropismN145_a":n145, "tropismN2_a":n2, "tropismN3_a":n3,
                "tropismS145_a":s145, "tropismS2_a":s2, "tropismS3_a":s3,
                "theta2_a": theta2, "theta3_a": theta3
                }
        if "_free" in exp_name:
            mods["bot_bc"] = "freeDrainage"  # otherwise potential
        if "_200" in exp_name:
            mods["water_table"] = 200  # otherwise 120
        run_soybean(exp_name, enviro_type, sim_time, mods, exp_name[:-4] + "/", save_all = True)

    elif type == "radii":
        print("running radii", flush = True)
        a145, a2, a3 = float(argv[5]), float(argv[6]), float(argv[7])
        hairsZone145, hairsZone2, hairsZone3 = float(argv[8]), float(argv[9]), float(argv[10])
        hairsLength145, hairsLength2, hairsLength3 = float(argv[11]), float(argv[12]), float(argv[13])
        mods = {
            "conductivity_mode": "scale",
            "scale_kr":1., "scale_kx":1.,
            "a145_a":a145, "a2_a":a2, "a3_a":a3,
            "hairsZone145_a":hairsZone145, "hairsZone2_a":hairsZone145, "hairsZone3_a":hairsZone145,
            "hairsLength145_a":hairsLength145, "hairsLength2_a":hairsLength2, "hairsLength3_a":hairsLength3
            }
        if "_free" in exp_name:
            mods["bot_bc"] = "freeDrainage"  # otherwise potential
        if "_200" in exp_name:
            mods["water_table"] = 200  # otherwise 120
        run_soybean(exp_name, enviro_type, sim_time, mods, exp_name[:-4] + "/", save_all = True)

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

        run_soybean(exp_name, enviro_type, sim_time, params, exp_name[:-4] + "/", save_all = True)

    else:

        print("Unknown run sa type", type)


if __name__ == "__main__":

    # enviro_type = 0
    # sim_time = 1
    # mods = {
    #     "output_times": [40],
    #     "conductivity_mode": "scale",
    #     "scale_kr":1.,
    #     "scale_kx":1.,
    #     }
    # run_soybean("test", enviro_type, sim_time, mods, save_all = True)
    # dd

    run(sys.argv)

    # elif type == "conductivities10":  # makes little sense to treat kr, kx independently
    #     print("running conductivities10", flush = True)
    #     mods = {
    #         "conductivity_mode": "age_dependent",
    #         "kr_young1": float(argv[5]),
    #         "kr_old1": float(argv[6]),
    #         "kr_young2": float(argv[7]),
    #         "kr_old2": float(argv[8]),
    #         "kr3": float(argv[9]),
    #         "kx_young1": float(argv[10]),
    #         "kx_old1": float(argv[11]),
    #         "kx_young2": float(argv[12]),
    #         "kx_old2": float(argv[13]),
    #         "kx3": float(argv[14])
    #     }
    #     run_soybean(exp_name, enviro_type, sim_time, mods, save_all = True)
    #
    # elif type == "singleroot_conductivities10":  # a single root example
    #     print("singleroot_conductivities10", flush = True)
    #     mods = {
    #         "conductivity_mode": "age_dependent",
    #         "kr_young1": float(argv[5]),
    #         "kr_old1": float(argv[6]),
    #         "kr_young2": float(argv[5]),
    #         "kr_old2": float(argv[6]),
    #         "kr3": float(argv[5]),
    #         "kx_young1": float(argv[7]),
    #         "kx_old1": float(argv[8]),
    #         "kx_young2": float(argv[7]),
    #         "kx_old2": float(argv[8]),
    #         "kx3": float(argv[7]),
    #         "filename": "data/Glycine_max_Moraes2020_singleroot.xml",
    #         "initial_age": 1.,
    #         "initial_totalpotential":-1000,
    #         "domain_size": [1., 1., 200.],
    #         "bot_bc": "noFlux",
    #     }
    #     run_soybean(exp_name, enviro_type, sim_time, mods, save_all = True)
