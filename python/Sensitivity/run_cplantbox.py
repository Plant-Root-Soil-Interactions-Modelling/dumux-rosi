"""
    starts a single root architecture simulation of soybean (indlucing Krs, SUF), 
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

import copy
import json
import timeit
import numpy as np
import matplotlib.pyplot as plt

import plantbox as pb
import scenario_setup as scenario
import hydraulic_model


class NpEncoder(json.JSONEncoder):

    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.bool):
            return bool(obj)
        return super(NpEncoder, self).default(obj)


def run_soybean(exp_name, enviro_type, sim_time, mods, save_all = False):
    """
        exp_name                experiment name (name for output files)
        enviro_type             envirotype (number)
        sim_time                simulation time [days]
    """

    if not (("conductivity_mode" in mods) and ("output_times" in mods)):
            print("run_cplantbox.run_soybean() conductivity_mode and output_times must be defined in mods dictionary")
            raise()

    sim_params = {"exp_name": exp_name, "enviro_type": enviro_type, "sim_time":sim_time}

    print("***********************")
    print("run_cplantbox.run_soybean", exp_name, enviro_type, sim_time)
    print(mods)
    print("***********************\n", flush = True)
    mods_copy = copy.deepcopy(mods)  # for json dump
    mods_copy.update(sim_params)

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

    r, params = hydraulic_model.create_mapped_rootsystem(min_b, max_b, cell_number, None, xml_name, stochastic = False, mods = mods, model = "Meunier")

    cm = mods["conductivity_mode"]
    if cm == "age_dependent":
        scenario.init_lupine_conductivities_sa(r.params, mods["kr_young1"], mods["kr_old1"], mods["kr_young2"], mods["kr_old2"], mods["kr3"],
                                               mods["kx_young1"], mods["kx_old1"], mods["kx_young2"], mods["kx_old2"], mods["kx3"])
        mods.pop("conductivity_mode")
        mods.pop("kr_young1")
        mods.pop("kr_old1")
        mods.pop("kr_young2")
        mods.pop("kr_old2")
        mods.pop("kr3")
        mods.pop("kx_young1")
        mods.pop("kx_old1")
        mods.pop("kx_young2")
        mods.pop("kx_old2")
        mods.pop("kx3")
    elif cm == "scale":
        scenario.init_lupine_conductivities(r.params, mods["scale_kr"], mods["scale_kx"])
        mods.pop("conductivity_mode")
        mods.pop("scale_kr")
        mods.pop("scale_kx")
    else:
        raise "run_cplantbox.run_soybean() conductivity_mode unknown"

    output_times = mods["output_times"]
    mods.pop("output_times")

    if mods:  # something unused in mods
        print("\n********************************************************************************")
        print("scenario_setup.create_mapped_rootsystem() WARNING unused parameters:")
        for k, v in mods.items():
            print("key:", k)
        print("********************************************************************************\n")
        # raise

    length, surface, volume, depth, RLDmean, RLDz, krs, SUFz, RLD, SUF = run_(r, sim_time, initial_age, area, output_times)

    print("writing", exp_name)
    # results
    scenario.write_cplantbox_results(exp_name, length, surface, volume, depth, RLDmean, RLDz, krs, SUFz, RLD, SUF, write_all = save_all)
    # parameters
    with open("results_cplantbox/" + exp_name + "_mods.json", "w+") as f:
        json.dump(mods_copy, f, cls = NpEncoder)
    # resulting cplantbox parameters
    # r.ms.writeParameters("results_cplantbox/" + exp_name + ".xml")  # corresponding xml parameter set

    print("finished " + exp_name)

    return r


def run_(r, sim_time, initial_age, area, out_times):
    """ starts a cplantbox simulation without dynamic soil """

    start_time = timeit.default_timer()

    out = out_times.copy()
    out.insert(0, 0.)
    out.append(sim_time)  # always export final simulation time
    dt_ = np.diff(out)
    length, surface, volume, depth, RLDmean, RLDz, krs, SUFz, RLD, SUF = [], [], [], [], [], [], [], [], [], []

    print("dt_", dt_)
    t = initial_age
    for dt in dt_:

        r.ms.simulate(dt, False)
        t += dt  # current simulation time
        length.append(np.sum(r.ms.getParameter("length")))
        volume.append(np.sum(r.ms.getParameter("volume")))  # we ignore root hairs
        surface.append(np.sum(r.ms.getParameter("surface")))  # TODO: include root hairs
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

    print ("\nrun_cplantbox.run_(): Root architecture simulation solved in ", timeit.default_timer() - start_time, " s")

    return np.array(length), np.array(surface), np.array(volume), np.array(depth), np.array(RLDmean), np.array(RLDz), np.array(krs), np.array(SUFz), np.array(RLD), np.array(SUF)


if __name__ == "__main__":

    enviro_type = 0
    sim_time = 87.5
    theta1 = None
    src = None
    kx = [0.1, 1.e-3, 1.e-3]  # cm3/day
    kx_old = [0.35, 0.015]
    kr = [1.e-3, 4.e-3, 4.e-3]  # 1/day
    kr_old = [5e-4, 0.0015]

    mods = {"output_times": [40], "conductivity_mode": "age_dependent",
            "kr_young1": kr[0], "kr_old1": kr_old[0], "kr_young2": kr[1], "kr_old2": kr_old[1], "kr3": kr[2],
            "kx_young1": kx[0], "kx_old1": kx_old[0], "kx_young2": kx[1], "kx_old2": kx_old[1], "kx3": kx[2]}
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
