"""
(experimental) plots data from global optimization 
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

import os
import json
import numpy as np
import matplotlib.pyplot as plt

import global_optimization_tools as got


def get_(name, all):
    l = np.array([a[name] for a in all])
    return l


def plot_1D(nameX, nameY, all):
    x = get_(nameX, all)
    y = get_(nameY, all)
    I = np.argsort(x)
    plt.plot(x[I], y[I])
    plt.xlabel(nameX)
    plt.ylabel(nameY)
    plt.show()


def scatter_1D(nameX, nameY, all, nameC = None):
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


""" font sizes """
SMALL_SIZE = 12
MEDIUM_SIZE = 16
BIGGER_SIZE = 16
plt.rc('font', size = SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize = SMALL_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize = MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize = SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize = BIGGER_SIZE)  # fontsize of the figure title

""" """
folder_path = "results_cplantbox/"
exp_name = "soybean_all14"

""" load everything & merge npz results into input parameter json"""

all = got.load_json_files(exp_name, folder_path)
got.merge_results(folder_path, all)

# scatter_1D("lmax1_a", "length", all)
# scatter_1D(["r1_a"], ["length"], all)
scatter_1D(["r145_a", "ln145_a"], ["length", "volume"], all)

# scatter_1D("hairsLength_a", "krs", all,  "lmax1_a")
# scatter_1D("hairsZone_a", "krs", all,  "lmax1_a")

