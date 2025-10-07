"""
    Macroscopic:

    re-runs an experiment from a json file as a root architecture simulation of soybean using run_cplantbox,
    if params["conductivity_mode"] == "from_mecha", set the mecha result path in L29
    
    Daniel Leitner, 2025       
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

import os
import json
import numpy as np

import plantbox as pb
import visualisation.vtk_plot as vp

import run_cplantbox
import scenario_setup
import carbon_cost


def show_me_file(exp_name, folder_path, vis_hairs = True, mode = "interactive", png_name = None):
    """ opens a json file and uses show_me() to run and plot simulation """

    """ open parameters """
    file_path = os.path.join(folder_path, exp_name + "_mods.json")
    with open(file_path, 'r', encoding = 'utf-8') as file:
        params = json.load(file)
    params["mecha_path"] = "/home/daniel/Dropbox/Code/granar/mecha_results"  # path changed

    """ open corresponding macroscopic results"""
    results = {}
    filename = params["exp_name"]
    assert exp_name == filename, "show_me(): there is some mix up with the filenames"
    filename += ".npz"
    file_path = os.path.join(folder_path, filename)
    results_ = np.load(file_path)  #
    i = -1  # last time point
    results["length"] = results_["length"][i]
    results["surface"] = results_["surface"][i]
    results["volume"] = results_["volume"][i]
    results["depth"] = np.abs(results_["depth"][i])  # they are stored with a (wrong) sign
    results["RLDmean"] = results_["RLDmean"][i]
    results["RLDz"] = np.abs(results_["RLDz"][i])  # they are stored with a (wrong) sign
    results["krs"] = results_["krs"][i]
    results["SUFz"] = np.abs(results_["SUFz"][i])  # they are stored with a (wrong) sign
    results["area"] = np.abs(results_["area"][i])
    print(results)

    show_me(params, mode, vis_hairs, png_name)


def show_me(params, mode = "interactive", vis_hairs = True, png_name = None):
    """ runs and plots simulation """

    """ rerun simulations (do not overwrite results) """
    enviro_type = params["enviro_type"]
    sim_time = params["sim_time"]
    params.pop("exp_name")
    params.pop("enviro_type")
    params.pop("sim_time")

    rs = run_cplantbox.run_soybean("dummy", enviro_type, sim_time, params, False)
    print("Carbon cost")
    print("volume    ", carbon_cost.carbon_cost(rs.ms, params, "volume"), "g")
    print("simple    ", carbon_cost.carbon_cost(rs.ms, params, "simple"), "g")
    print("anatomical", carbon_cost.carbon_cost(rs.ms, params, "anatomical"), "g")

    # root hair visualization hack
    ms = rs.ms
    if vis_hairs:
        n = len(ms.radii)
        radii = np.array([ms.getEffectiveRadius(i) for i in range(0, n)])
        ms.radii = radii

    ana = pb.SegmentAnalyser(ms.mappedSegments())

    if mode == "interactive":
        vp.plot_plant(ana, "subType")
    elif mode == "png":
        actors, color_bar = vp.plot_plant(ana, "age", render = False)
        tube_plot_actor = actors[0]
        actor = actors[1]
        ren = vp.render_window([tube_plot_actor, actor], "plot_plant", color_bar, tube_plot_actor.GetBounds(), True)
        if png_name == None:
            png_name = exp_name
        vp.write_png(ren.GetRenderWindow(), png_name)
    else:
        raise "show_me(): unknown mode"


def compare(exp1_name, exp2_name, folder_path):
    """ open parameters """
    file_path = os.path.join(folder_path, exp1_name + "_mods.json")
    with open(file_path, 'r', encoding = 'utf-8') as file:
        params1 = json.load(file)
    params1["mecha_path"] = "/home/daniel/Dropbox/Code/granar/mecha_results"  # path changed
    file_path = os.path.join(folder_path, exp2_name + "_mods.json")
    with open(file_path, 'r', encoding = 'utf-8') as file:
        params2 = json.load(file)
    params2["mecha_path"] = "/home/daniel/Dropbox/Code/granar/mecha_results"  # path changed

    scenario_setup.attach_conductivitis(params1)
    scenario_setup.attach_conductivitis(params2)

    print()
    keys_of_interest = ["lmax145_a", "r145_a", "ln145_a", "a145_a"]
    for k in keys_of_interest:
        print(k, params1[k], params2[k], np.abs((params1[k] - params2[k]) / max(params1[k], params2[k])))

    print()
    keys_of_interest = ["lmax2_a", "r2_a", "ln2_a", "a2_a" ]
    for k in keys_of_interest:
        print(k, params1[k], params2[k], np.abs((params1[k] - params2[k]) / max(params1[k], params2[k])))

    print()
    keys_of_interest = ["lmax3_a", "r3_a", "a3_a"]
    for k in keys_of_interest:
        print(k, params1[k], params2[k], np.abs((params1[k] - params2[k]) / max(params1[k], params2[k])))

    filename1 = params1["exp_name"]
    assert exp1_name == filename1, "show_me(): there is some mix up with the filenames"
    filename1 += ".npz"
    file_path = os.path.join(folder_path, filename1)
    results1_ = np.load(file_path)
    filename2 = params2["exp_name"]
    assert exp2_name == filename2, "show_me(): there is some mix up with the filenames"
    filename2 += ".npz"
    file_path = os.path.join(folder_path, filename2)
    results2_ = np.load(file_path)
    print()
    print("depth", np.abs(results1_["depth"][-1]), np.abs(results2_["depth"][-1]))
    print("volume", np.abs(results1_["volume"][-1]), np.abs(results2_["volume"][-1]))


if __name__ == "__main__":

    folder_path = "results_cplantbox/"
    exp_name = "soybean_all14_64b677a63373c8db3267054d44828238b00525ebb83abca8c85930e130efaaa6"
    exp_name1 = "soybean_all14_3564a0f636e9f600fe68bf96ffca4124135106ae4787b9a3332334b04abcdf1a"  # Number 213 from cluster (1,0) -> see analyze_transpiration_cluster.py
    exp_name2 = "soybean_all14_f7319dc5d83c72932fd39e4afbf6e50822c2f7bf13b27fc5749c1128642a95d2"  # Number 138 from cluster (1,0) -> see analyze_transpiration_cluster.py
    show_me_file(exp_name2, folder_path)
    # compare(exp_name1, exp_name2, folder_path)

