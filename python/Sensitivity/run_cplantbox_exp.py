"""
    starts a single root architecture simulation of soybean (indlucing Krs, SUF), 
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

import os
import json
import numpy as np

import visualisation.vtk_plot as vp

import run_cplantbox


def show_me(exp_name, folder_path):
    """ open parameter modifications """
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
    print(results)

    """ rerun simulations (do not overwrite results) """
    enviro_type = params["enviro_type"]
    sim_time = params["sim_time"]
    params.pop("exp_name")
    params.pop("enviro_type")
    params.pop("sim_time")
    rs = run_cplantbox.run_soybean("dummy", enviro_type, sim_time, params, False)

    vp.plot_plant(rs.ms, "age")


if __name__ == "__main__":

    folder_path = "results_cplantbox/"
    exp_name = "soybean_length14_6a8fba2616c035a3091f02f27a172cfdda8f89861027dcf177edc55776e81028"

    show_me(exp_name, folder_path)
