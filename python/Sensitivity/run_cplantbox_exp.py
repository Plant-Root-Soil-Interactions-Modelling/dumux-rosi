"""
    reruns an experiment from a json file as a root architecture simulation of soybean (indlucing Krs, SUF), 
    using run_cplantbox
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

import os
import json
import numpy as np

import plantbox as pb
import visualisation.vtk_plot as vp

import run_cplantbox


def show_me_file(exp_name, folder_path, vis_hairs = True, mode = "interactive", png_name = None):
    """ opens a json file and uses show_me() to run and plot simulation """

    """ open parameters """
    file_path = os.path.join(folder_path, exp_name + "_mods.json")
    with open(file_path, 'r', encoding = 'utf-8') as file:
        params = json.load(file)
    print(params)

    """ open corresponding results """
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
    """ runs and plots simulaiton """

    """ rerun simulations (do not overwrite results) """
    enviro_type = params["enviro_type"]
    sim_time = params["sim_time"]
    params.pop("exp_name")
    params.pop("enviro_type")
    params.pop("sim_time")

    rs = run_cplantbox.run_soybean("dummy", enviro_type, sim_time, params, False)

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
        actors, color_bar = vp.plot_plant(ana, "subType", render = False)
        tube_plot_actor = actors[0]
        actor = actors[1]
        ren = vp.render_window([tube_plot_actor, actor], "plot_plant", color_bar, tube_plot_actor.GetBounds(), True)
        if png_name == None:
            png_name = exp_name
        vp.write_png(ren.GetRenderWindow(), png_name)
    else:
        raise "show_me(): unknown mode"


if __name__ == "__main__":

    folder_path = "results_cplantbox/"
    exp_name = "soybean_all14_64b677a63373c8db3267054d44828238b00525ebb83abca8c85930e130efaaa6"

    show_me_file(exp_name, folder_path)
