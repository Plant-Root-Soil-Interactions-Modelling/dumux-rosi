"""
(experimental) plots data from global optimization 
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

import os
import json
import numpy as np
import matplotlib.pyplot as plt
# import pandas as pd


def load_json_files(exp_name, folder_path):
    """ opens all json files starting with exp_name in a folder and returns a list of Python dictionaries 
        (one dictionary by file)"""
    json_data_list = []
    for filename in os.listdir(folder_path):
        if filename.startswith(exp_name) and filename.endswith(".json"):
            file_path = os.path.join(folder_path, filename)
            try:
                with open(file_path, 'r', encoding = 'utf-8') as file:
                    data = json.load(file)
                    json_data_list.append(data)
            except (json.JSONDecodeError, IOError) as e:
                print(f"Error loading {filename}: {e}")

    print("global_optimization_tools.load_json_files(): opened", len(json_data_list), "files")
    return json_data_list


def merge_results(folder_path, json_data_list, i = -1):
    """ Merges the json_data (list of dict) with simulation targets given within .npz with corresponding name """
    for data in json_data_list:

        filename = data["exp_name"]
        filename += ".npz"
        file_path = os.path.join(folder_path, filename)
        try:
            results = np.load(file_path)
            data["length"] = results["length"][i]
            data["surface"] = results["surface"][i]
            data["volume"] = results["volume"][i]
            data["depth"] = np.abs(results["depth"][i])  # they are stored with a (wrong) sign
            data["RLDmean"] = results["RLDmean"][i]
            data["RLDz"] = np.abs(results["RLDz"][i])  # they are stored with a (wrong) sign
            data["krs"] = results["krs"][i]
            data["SUFz"] = np.abs(results["SUFz"][i])  # they are stored with a (wrong) sign

        except (json.JSONDecodeError, IOError) as e:
            print(f"Error loading {filename}: {e}")


def get_(name, all):
    """ returns the key @param name from a list of dicts as np.array """
    l = np.array([a[name] for a in all])
    return l


def fetch_features(feature_names, all):
    """ returns feature_names as np.array with the shape (len(data), len(features)) """
    data = []
    for name in feature_names:
        data.append(get_(name, all))

    data = np.array(data).transpose()
    return data


def print_info(data, feature_names):
    """ prints features min, max, mean, and sd (featrure names labels the data columns) """
    print("name\tmin\t\tmax\t\tmean\t\tstd")
    for i, name in enumerate(feature_names):
        # print(name, "\t", np.min(data[:, i]), "\t", np.max(data[:, i]), "\t", np.mean(data[:, i]), "\t", np.std(data[:, i]))
        print("{:s}\t{:10.4f}\t{:10.4f}\t{:10.4f}\t{:10.4f}".format(name, np.min(data[:, i]), np.max(data[:, i]), np.mean(data[:, i]), np.std(data[:, i])))


def filter_data(data, feature_ind, min_, max_):
    """ crops data to where feature_ind stays within (min_, max_) """
    ind_ = np.bitwise_and(data[:, feature_ind] > min_, data[:, feature_ind] < max_)
    data = data[ind_,:]
    return data


def scale_data(data):
    """ scales data to have a mean of 0 and a std of 1 """
    mean = np.mean(data, axis = 0)
    std = np.std(data, axis = 0)
    data = data - mean
    data = np.divide(data, std)
    return data


def plot_1D(nameX, nameY, all):
    """ plots nameX vs nameY sorting the values along the x-axis """
    x = get_(nameX, all)
    y = get_(nameY, all)
    I = np.argsort(x)
    plt.plot(x[I], y[I])
    plt.xlabel(nameX)
    plt.ylabel(nameY)
    plt.show()


def scatter_1D(nameX, nameY, all, nameC = None):
    """ scatter plots one or multiple pairs of nameX vs nameY """
    if isinstance(nameX, str):
        nameX = [nameX]
    if isinstance(nameY, str):
        nameY = [nameY]
    assert len(nameX) == len(nameY), "oh no"

    n = np.ceil(np.sqrt(len(nameX)))
    print(n, np.sqrt(len(nameX)))
    fig, ax = plt.subplots(int(n), int(n), figsize = (16, 16))
    if n == 1:
        ax = np.array([ax])
    ax_ = ax.flat

    for i, _ in enumerate(nameX):
        x = get_(nameX[i], all)
        y = get_(nameY[i], all)
        if nameC is not None:
            c = get_(nameC, all)
        I = np.argsort(x)
        if nameC is not None:
            ax_[i].scatter(x[I], y[I], c = c[I])
        else:
            ax_[i].scatter(x[I], y[I])
        ax_[i].set_xlabel(nameX[i])
        ax_[i].set_ylabel(nameY[i])

    plt.show()


if __name__ == "__main__":

    folder_path = "results_cplantbox/"
    exp_name = "soybean_test"

    all = load_json_files(exp_name, folder_path)
    merge_results(folder_path, all)

    # scatter_1D("lmax1_a", "length", all)
    # scatter_1D(["r1_a"], ["length"], all)
    # scatter_1D(["r1_a", "ln1_a"], ["length", "volume"], all)

    # scatter_1D("hairsLength_a", "krs", all,  "lmax1_a")
    # scatter_1D("hairsZone_a", "krs", all,  "lmax1_a")

    print("data set", len(all), "entries")
    print(all[10])

    # StandardScaler (sklearn)

    # step 1 (scale features, not targets)
    # step 2 (pca, pca.explained_variance ratio

    feature_names = ['a1_a', 'lmax1_a', 'ln1_a', 'r1_a', 'tropismN1_a', 'hairsLength_a', 'hairsZone_a', 'hairsElongation_a']
    data = fetch_features(feature_names, all)
    data = scale_data(data)

