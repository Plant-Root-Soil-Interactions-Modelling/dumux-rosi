"""
    Macroscopic:
    
    Starts a single root architecture simulation of soybean (calculating Krs, and SUF),    
    
    see also run_sra.py (for the full dynamic simulation)
    
    Daniel Leitner, 2025        
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

import os
import copy
import json
import timeit
import numpy as np
import matplotlib.pyplot as plt

import plantbox as pb
import scenario_setup as scenario
import hydraulic_model
import carbon_cost


class NpEncoder(json.JSONEncoder):

    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.bool):
            return bool(obj)
        return super(NpEncoder, self).default(obj)


def run_soybean(exp_name, enviro_type, sim_time, params, save_all = False):
    """
        exp_name                experiment name (name for output files)
        enviro_type             envirotype (number)
        sim_time                simulation time [days]
        params                  parameters that are adjusted from the base parametrization (default: data/Glycine_max_Moraes2020_opt2_modified.xml)
        save_all                ---
    """

    mods = params.copy()

    if not (("conductivity_mode" in mods) and ("output_times" in mods)):
            print("run_cplantbox.run_soybean() conductivity_mode and output_times must be defined in mods dictionary")
            raise()

    sim_params = {"exp_name": exp_name, "enviro_type": enviro_type, "sim_time":sim_time}

    print("***********************")
    print("run_cplantbox.run_soybean", exp_name, enviro_type, sim_time)
    print(mods)
    print("***********************\n", flush = True)
    mods_copy = copy.deepcopy(mods)  # for json dump
    mods_copy.update(sim_params)  # all input arguments, to make reproducable

    if "filename" in mods:
        xml_name = mods["filename"]
        mods.pop("filename")
    else:
        xml_name = "data/Glycine_max_Moraes2020_opt2_modified.xml"  # root growth model parameter file

    if "initial_age" in mods:  # used also in hydraulic_model
        initial_age = mods["initial_age"]
    else:
        initial_age = 1.  # day

    soil_, table_name, min_b, max_b, cell_number, area, Kc = scenario.soybean(int(enviro_type))

    if "domain_size" in mods:
        domain_size = np.array(mods["domain_size"])
        min_b = -np.array([domain_size[0] / 2., domain_size[1] / 2., domain_size[2]])
        max_b = np.array([domain_size[0] / 2., domain_size[1] / 2., 0.])
        area = domain_size[0] * domain_size[1]
        # print("Changed domain size to", min_b, max_b, area)
        mods.pop("domain_size")

    cdata = scenario.prepare_conductivities(mods)
    r, hydraulic_params = hydraulic_model.create_mapped_rootsystem(min_b, max_b, cell_number, None, xml_name, stochastic = False, mods = mods, model = "Meunier")
    scenario.set_conductivities(hydraulic_params, mods, cdata)

    output_times = mods["output_times"]
    mods.pop("output_times")

    if mods:  # something unused in mods
        print("\n********************************************************************************")
        print("scenario_setup.create_mapped_rootsystem() WARNING unused parameters:")
        for k, v in mods.items():
            print("key:", k)
        print("********************************************************************************\n")
        # raise

    length, surface, volume, depth, RLDmean, RLDz, krs, SUFz, RLD, SUF, area_, carbon = run_(r, sim_time, initial_age, area, output_times, params)

    # results
    print("writing", exp_name)
    scenario.write_cplantbox_results(exp_name, length, surface, volume, depth, RLDmean, RLDz, krs, SUFz, RLD, SUF, area_, carbon, write_all = save_all)
    with open("results_cplantbox/" + exp_name + "_mods.json", "w+") as f:
        json.dump(mods_copy, f, cls = NpEncoder)
    # r.ms.writeParameters("results_cplantbox/" + exp_name + ".xml")  # corresponding xml parameter set

    print("finished " + exp_name)

    return r


def run_(r, sim_time, initial_age, area, out_times, params):
    """ starts a cplantbox simulation without dynamic soil """
    start_time = timeit.default_timer()

    out = out_times.copy()
    out.insert(0, 0.)
    out.append(sim_time)  # always export final simulation time
    dt_ = np.diff(out)
    length, surface, volume, depth, area_, RLDmean, RLDz, krs, SUFz, RLD, SUF, carbon = [], [], [], [], [], [], [], [], [], [], [], []

    print("dt_", dt_)
    t = initial_age
    for dt in dt_:

        r.ms.simulate(dt, False)
        t += dt  # current simulation time
        length.append(np.sum(r.ms.getParameter("length")))
        volume.append(np.sum(r.ms.getParameter("volume")))  # we ignore root hairs
        suf = r.get_suf(t)
        krs_, _ = r.get_krs(t)
        krs.append(krs_)
        z = np.linspace(0., -200, 200)
        sa = pb.SegmentAnalyser(r.ms.mappedSegments())
        sa.addData("suf", suf)
        rld = np.array(sa.distribution("length", 0., -200., 200, False)) / (area * 1)  # divide through layer volume
        RLD.append(rld)
        rld_ = rld / np.sum(rld)  # normalize
        RLDz.append(rld_.dot(z))
        RLDmean.append(np.mean(RLD))
        suf = np.array(sa.distribution("suf", 0., -200., 200, False))
        SUF.append(suf)
        SUFz.append(suf.dot(z))
        depth.append(sa.getMinBounds().z)
        bounds = sa.getMaxBounds().minus(sa.getMinBounds())
        area_.append(bounds.x * bounds.y)
        n = len(r.ms.radii)
        radii = np.array([r.ms.getEffectiveRadius(i) for i in range(0, n)])  # length including root hairs
        lengths = np.array(r.ms.segLength())
        surface.append(np.sum(2.*np.pi * np.multiply(radii, lengths)))
        c = carbon_cost.carbon_cost(r.ms, params, model = "simple")
        c = carbon_cost.carbon_cost(r.ms, params, model = "anatomical")
        carbon.append(c)

    print ("\nrun_cplantbox.run_(): Root architecture simulation solved in ", timeit.default_timer() - start_time, " s")

    return np.array(length), np.array(surface), np.array(volume), np.array(depth), np.array(RLDmean), np.array(RLDz), np.array(krs), np.array(SUFz), np.array(RLD), np.array(SUF), np.array(area_), np.array(carbon)


if __name__ == "__main__":

    enviro_type = 0
    sim_time = 87.5
    theta1 = None
    src = None

    mods = {"output_times": [40],
            "conductivity_mode": "from_mecha",
            "mecha_path": "/home/d.leitner/Dropbox/Code/granar/mecha_results",
            #"mecha_path": "~/Dropbox/Code/granar/mecha_results",
            'conductivity_index1': 1,
            'conductivity_index2': 30,
            'conductivity_index3': 100,
            "conductivity_age1": 14,
            "conductivity_age2": 7,
            "conductivity_age3": 7,
            }
    run_soybean("test", enviro_type, sim_time, mods, save_all = True)

    # type = sys.argv[1]
    # exp_name = sys.argv[2]
    # enviro_type = int(float(sys.argv[3]))
    # sim_time = float(sys.argv[4])
    #
    # if type == "original":
    #
    #     print("running original sa", flush = True)
    #
    #     kr = float(sys.argv[5])
    #     kx = float(sys.argv[6])
    #     lmax1 = float(sys.argv[7])
    #     lmax2 = float(sys.argv[8])
    #     lmax3 = float(sys.argv[9])
    #     theta1 = float(sys.argv[10])
    #     r1 = float(sys.argv[11])
    #     r2 = float(sys.argv[12])
    #     a = float(sys.argv[13])
    #     src = int(float(sys.argv[14]))
    #
    #     mods = {"lmax145":lmax1, "lmax2":lmax2, "lmax3":lmax3, "r145":r1, "r2":r2, "r3":r2, "a":a}
    #     mods["theta45"] = theta1
    #     mods["src"] = src
    #
    #     run_soybean(exp_name, enviro_type, sim_time, mods, kr, kx, save_all = True)
