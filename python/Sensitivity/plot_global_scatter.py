"""
    Macroscopic:
    
    scatter plots data from global static hydraulic model optimization, script form L188-
    
    Daniel Leitner, 2025    
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

import os
import zipfile
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

max_ind = 1
pbounds = {
    'conductivity_index1': (0, max_ind),
    'conductivity_index2': (0, max_ind),
    'conductivity_index3': (0, max_ind),
    "conductivity_age1": (1, 21),
    "conductivity_age2": (1, 21),
    "conductivity_age3": (1, 21),
    'src_a': (3, 11),
    'src_first_a': (3, 14),
    'src_delay_a': (3, 14),
    'lmax145_a': (50, 150),
    'ln145_a': (0.5, 10.),
    'r145_a': (0.2, 7.),
    'theta145_a': (np.pi / 8., np.pi / 2.),
    'tropismN145_a': (0., 3.5),
    'hairsLength145_a': (0., 0.1),
    'hairsZone145_a': (0., 5.),
    'lmax2_a': (5., 50.),
    'ln2_a': (0.5, 10.),
    'r2_a': (0.2, 7.),
    'theta2_a': (np.pi / 8., np.pi / 2.),
    'tropismN2_a': (0., 3.5),
    'hairsLength2_a': (0., 0.1),
    'hairsZone2_a': (0., 10.),
    'lmax3_a': (0.5, 50.),
    'r3_a': (0.2, 7.),
    'theta3_a': (np.pi / 8., np.pi / 2.),
    'tropismN3_a': (0., 3.5),
    'hairsLength3_a': (0., 0.1),
    'hairsZone3_a': (0., 10.),
    }
pbounds_depth = {
    'src_a': (3, 11),
    'src_first_a': (3, 14),
    'src_delay_a': (3, 14),
    'lmax145_a': (50, 150),
    'ln145_a': (0.5, 10.),
    'r145_a': (0.2, 7.),
    'theta145_a': (np.pi / 8., np.pi / 2.),
    'tropismN145_a': (0., 3.5),
    'lmax2_a': (5., 50.),
    'ln2_a': (0.5, 10.),
    'r2_a': (0.2, 7.),
    'theta2_a': (np.pi / 8., np.pi / 2.),
    'tropismN2_a': (0., 3.5),
    'lmax3_a': (0.5, 50.),
    'r3_a': (0.2, 7.),
    'theta3_a': (np.pi / 8., np.pi / 2.),
    'tropismN3_a': (0., 3.5),
    "anatomy_index145":0,
    "anatomy_index2":0,
    "anatomy_index3":0,
    # 'kr_young145': (0.2, 7.),
    # 'kr_young2': (np.pi / 8., np.pi / 2.),
    # 'kr_young3': (0., 3.5),
    # "kx_old145":0,
    # "kx_old2":0,
    # "kx_old3":0,
    }

""" """
folder_path = "results_cplantbox/"
# exp_name = "soybean_all14"
# exp_name = "soybean_some14"
# exp_name = "soybean_lmax14"
exp_name = "soybean_all14"

""" load everything & merge npz results into input parameter json"""
# all = got.load_json_files(exp_name, folder_path)
# got.merge_results(folder_path, all, 0)  # check signs for new simulations
# # print(all[0]["output_times"])

# Read JSON data from the ZIP file (faster)
with zipfile.ZipFile(exp_name + ".zip", "r") as zipf:
    with zipf.open(exp_name + ".json", "r") as json_file:
        all = json.load(json_file)  # Deserialize JSON data
all = list(all.values())

print("Number of simulations", len(all))
# print(all[0].keys())

all = got.filter_list(all, "length", 200., 30000)  # 76 * 3  *100 * 0.6 = 13680 cm;
print("Number after filtering", len(all))

""" restrict scatter to pareto (uncomment to use all) """
cluster_target_names = ["-carbon", "surface", "depth", "krs", "SUFz", "RLDz"]
ind = got.pareto_list(all, cluster_target_names)
pareto = []
for i, a in enumerate(all):
    all[i]["pareto"] = ind[i]
    if ind[i] == True:
        pareto.append(a)
all = pareto
print("Number of Pareto solutions", len(all))

""" perform clustering on scaled targets """
data = got.fetch_features(cluster_target_names, all)  # list of dicts with parameter sets and results
scaled_data = got.scale_data(data)
k = 9
centers, dist = kmeans(scaled_data, k, rng = 1)
# sample2node, dist = vq(scaled_data, centers)
node2sample, sample2node, node2sample_ = got.make_kmeans_maps(scaled_data, centers, n = 1)
for i, a in enumerate(all):
    a["node"] = sample2node[i]
    a["id"] = i

# part10m = []
# for i in all:  #
#     if i["length"] > 9000 and i["length"] < 11000:
#         part10m.append(i)
# print("part10m", len(part10m))

target_names = ["carbon", "krs", "SUFz", "surface", "RLDz"]  # "surface", "length", "volume", "depth", "RLDz", "krs", "SUFz", "RLDz", "SUFz"

# key_.extend(["kr1", "kx1", "ykr1", "ykx1", "a1", "kr2", "kx2", "ykr2", "ykx2", "a2", "kr3", "kx3", "ykr3", "ykx3", "a3"])

# got.scatter_1D(["lmax145_a", "ln145_a", "r145_a", "lmax2_a", "ln2_a", "r2_a", "lmax3_a", "r3_a", "src_a", "src_delay_a"],
#                 ["depth"] * 10, all, nameC = "id")  # , all, nameC = "id" # part10m

for target_name in  target_names:

    key_ = list(pbounds_depth.keys())  # [9:]  # skip the first three

    slopes = got.scatter_1D(key_, [target_name] * len(key_), all, nameC = "node", scale = True, show = False)
    I = np.argsort(-np.abs(slopes))

    print("\n")  # impact order
    print(target_name)
    for i in I[:5]:
        print(key_[i], f"{slopes[i]:g}")

    plt.savefig("figures/scatter_" + target_name + ".png", dpi = 150, bbox_inches = 'tight')
    print("saved: " + " figures/scatter_" + target_name + ".png")

# key_ = ["kr1", "kx1", "ykr1", "ykx1", "a1", "kr2", "kx2", "ykr2", "ykx2", "a2", "kr3", "kx3", "ykr3", "ykx3", "a3"]
# slopes = got.scatter_1D(key_, target_names * len(key_), all, nameC = "id", scale = True)  # , all, nameC = "id" # part10m
# I = np.argsort(-np.abs(slopes))
# print("\n")  # impact order
# for i in I:
#     print(key_[i], slopes[i])

# got.scatter_1D(["lmax145_a", "ln145_a", "src_a", "src_delay_a"], ["length"] * 4, all, nameC = "id")

# got.scatter_1D_cross(["lmax145_a", "ln145_a", "r145_a", "lmax2_a", "ln2_a", "r2_a", "lmax3_a", "r3_a", "src_a", "src_delay_a"], part10m)
