"""
    scatter plots data from global static hydraulic model optimization 
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

import os
import json
import numpy as np
from scipy.cluster.vq import vq, whiten, kmeans
import matplotlib.pyplot as plt

from minisom import MiniSom

import global_optimization_tools as got

""" font sizes """
SMALL_SIZE = 12
MEDIUM_SIZE = 16
BIGGER_SIZE = 16
plt.rc('font', size = SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize = SMALL_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize = MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize = BIGGER_SIZE)  # legend fontsize
plt.rc('figure', titlesize = BIGGER_SIZE)  # fontsize of the figure title


def distortion_plot(data):
    data_ = whiten(data)
    dist_ = []
    for i in range(0, 20):
        centers, dist = kmeans(data_, i + 1)
        dist_.append(dist)
    plt.plot(range(0, 20), dist_)
    plt.ylabel("distortion")
    plt.xlabel("number of clusters")
    plt.show()


""" """
folder_path = "results_cplantbox/"
# exp_name = "soybean_all14"
# exp_name = "soybean_some14"
# exp_name = "soybean_lmax14"
exp_name = "soybean_all14"

""" load everything & merge npz results into input parameter json"""
all = got.load_json_files(exp_name, folder_path)
got.merge_results(folder_path, all)  # check signs for new simulations
all = got.filter_list(all, "length", 100., 30000)  # 76 * 3  *100 * 0.6 = 13680 cm;

target_names = ["length"]  # "surface", "length", "volume", "depth", "RLDz", "krs", "SUFz", "RLDz", "SUFz"
data = got.fetch_features(target_names, all)  # list of dicts with parameter sets and results
print("\ndata")
print(data.shape)
scaled_data = got.scale_data(data)

""" perform clustering on scaled targets """
k = 9
centers, dist = kmeans(scaled_data, k, rng = 1)
# sample2node, dist = vq(scaled_data, centers)
node2sample, sample2node, node2sample_ = got.make_kmeans_maps(scaled_data, centers, n = 1)
for i, a in enumerate(all):
    a["node"] = sample2node[i]
    a["id"] = i
print("node of data 0")
print(all[0]["node"])

ind = node2sample_[4]  # mid
part4 = [all[i] for i in ind]

part10m = []
for i in all:  #
    if i["length"] > 9000 and i["length"] < 11000:
        part10m.append(i)
print("part10m", len(part10m))

got.scatter_1D(["lmax145_a", "ln145_a", "r145_a", "lmax2_a", "ln2_a", "r2_a", "lmax3_a", "r3_a", "src_a", "src_delay_a"],
                ["depth"] * 10, all, nameC = "id")  # , all, nameC = "id" # part10m

# got.scatter_1D(["lmax145_a", "ln145_a", "src_a", "src_delay_a"], ["length"] * 4, all, nameC = "id")

# got.scatter_1D_cross(["lmax145_a", "ln145_a", "r145_a", "lmax2_a", "ln2_a", "r2_a", "lmax3_a", "r3_a", "src_a", "src_delay_a"], part10m)
